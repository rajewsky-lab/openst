import argparse
import logging
import multiprocessing
import os
import shutil
import subprocess
import sys
import tempfile
import threading
import xml.etree.ElementTree as ET
from concurrent.futures import ProcessPoolExecutor
from typing import List, Tuple

from tqdm import tqdm

# for the last merge sort, in gigabytes
MEM_PER_CORE = 4

log_lock = threading.Lock()

def run_command(cmd: List[str], description: str) -> Tuple[int, str, str]:
    """Run a command and return its exit code, stdout, and stderr."""
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate()
    return process.returncode, stdout, stderr

def log_output(description: str, returncode: int, stdout: str, stderr: str):
    """Log the output of a command in a synchronized manner."""
    with log_lock:
        logging.debug(f"\n{'=' * 40}\n{description}\n{'=' * 40}")
        if stdout:
            logging.debug("STDOUT:\n%s", stdout)
        if stderr and returncode:
            logging.warning("STDERR:\n%s", stderr)
        logging.debug("Return code: %d\n", returncode)

def parse_run_info_xml(run_info_path: str) -> List[str]:
    """Parse RunInfo.xml and extract tile information."""
    tree = ET.parse(run_info_path)
    root = tree.getroot()

    tiles = []
    for tile_elem in root.findall(".//Tile"):
        tiles.append(tile_elem.text.strip())

    return tiles


def process_tile(tile: str, args: argparse.Namespace):
    bcl_out_dir = os.path.join(args.bcl_out, tile)
    fastq_file = os.path.join(bcl_out_dir, "Undetermined_S0_R1_001.fastq.gz")

    os.makedirs(bcl_out_dir, exist_ok=True)

    if args.demux_tool == "bcl2fastq":
        demux_cmd = ["bcl2fastq", "-R", args.bcl_in, "--no-lane-splitting", "-o", bcl_out_dir, "--tiles", f"s_{tile}"]
    elif args.demux_tool == "bcl-convert":
        demux_cmd = ["bcl-convert", "--bcl-input-directory", args.bcl_in, "--no-lane-splitting", "true", "--force", "--no-sample-sheet", "true", "--output-directory", bcl_out_dir, "--tiles", f"s_{tile}"]
    else:
        logging.error("The demultiplexing tool has to be 'bcl2fastq' or 'bcl-convert'")
        raise ValueError("The demultiplexing tool has to be 'bcl2fastq' or 'bcl-convert'")
    
    returncode, stdout, stderr = run_command(demux_cmd, f"{args.demux_tool} --tile {tile}")
    log_output(f"bcl2fastq for tile {tile}", returncode, stdout, stderr)
    if returncode != 0:
        logging.warning(f"return code {returncode}, at command '{demux_cmd}'")

    os.makedirs(args.tilecoords_out, exist_ok=True)

    openst_cmd = [
        "openst",
        "barcode_preprocessing",
        "--fastq-in",
        fastq_file,
        "--tilecoords-out",
        args.tilecoords_out,
        "--out-suffix",
        args.out_suffix,
        "--out-prefix",
        f"{args.out_prefix}{tile}",
        "--crop-seq",
        args.crop_seq,
        "--single-tile",
    ]
    if args.rev_comp:
        openst_cmd.append("--rev-comp")

    returncode, stdout, stderr = run_command(openst_cmd, f"openst barcode_preprocessing --tile {tile}")
    log_output(f"openst barcode_preprocessing for tile {tile}", returncode, stdout, stderr)
    if returncode != 0:
        logging.warning(f"return code {returncode}, at command '{openst_cmd}'")

def process_deduplication_file(file: str, input_folder: str, output_folder: str):
    input_file = os.path.join(input_folder, file)
    output_file = os.path.join(output_folder, f".uncompressed.sorted.{file}")
    cmd = f"""
    zcat "{input_file}" | 
    awk -F'\\t' 'NR==1{{print $0"\\tfilename"; next}} {{print $0"\\t{file}"}}' | 
    tail -n +2 | 
    sort -u -k1,1 -t$'\\t' > "{output_file}"
    """
    returncode, stdout, stderr = run_command(["bash", "-c", cmd], f"deduplicate --tile {file}")
    log_output(f"deduplicate --tile {file}", returncode, stdout, stderr)
    if returncode != 0:
        raise subprocess.CalledProcessError(returncode, cmd)

def deduplicate_tiles(input_folder: str, output_folder: str):
    os.makedirs(output_folder, exist_ok=True)
    files = [f for f in os.listdir(input_folder) if f.endswith('.txt.gz')]

    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_deduplication_file, file, input_folder, output_folder) for file in files]
        for _ in tqdm(futures, total=len(futures), desc="deduplicate"):
            _.result()

def merge_tiles(input_folder: str, output_file: str, threads: int):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    files = [os.path.join(input_folder, f) for f in os.listdir(input_folder) if f.startswith(".uncompressed.sorted.")]
    
    # Create a temporary file with the list of files to sort
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_file:
        for file in files:
            temp_file.write(f"{file}\n")
        temp_file_name = temp_file.name

    try:
        # Use process substitution to pass the file list to sort
        cmd = f"sort --parallel={threads} --batch-size=128 --buffer-size={max(4, int(MEM_PER_CORE*threads))}G -m -u -k1,1 -t $'\\t' -T {input_folder} $(cat {temp_file_name}) > {output_file}"
        returncode, stdout, stderr = run_command(["bash", "-c", cmd], "merge")
        log_output("merge", returncode, stdout, stderr)
        if returncode != 0:
            logging.warning(f"return code {returncode}, at command '{cmd}'")
    finally:
        # Clean up the temporary file
        os.unlink(temp_file_name)


def distribute_per_file(input_file: str, output_folder: str):
    os.makedirs(output_folder, exist_ok=True)
    cmd = f'awk -F\'\\t\' \'{{print $1"\\t"$2"\\t"$3 > "{output_folder}/"$4}}\' {input_file}'
    returncode, stdout, stderr = run_command(["bash", "-c", cmd], "distribute")
    log_output("distribute", returncode, stdout, stderr)
    if returncode != 0:
        logging.warning(f"return code {returncode}, at command '{cmd}'")


def process_compression_file(file: str, input_folder: str, output_folder: str):
    input_file = os.path.join(input_folder, file)
    output_file = os.path.join(output_folder, f"{file}.gz")
    cmd = f"gzip -c {input_file} > {output_file}"
    returncode, stdout, stderr = run_command(["bash", "-c", cmd], f"gzip {file}")
    log_output(f"gzip {file}", returncode, stdout, stderr)
    if returncode != 0:
        logging.warning(f"return code {returncode}, at command '{cmd}'")

def compress_files(input_folder: str, output_folder: str):
    os.makedirs(output_folder, exist_ok=True)
    files = [f for f in os.listdir(input_folder) if not f.startswith(".")]
    
    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_compression_file, file, input_folder, output_folder) for file in files]
        for _ in tqdm(futures, total=len(files), desc="gzip"):
            _.result()  # This will raise any exceptions that occurred during execution


def check_bcl2fastq():
    """Check if bcl2fastq is installed and accessible."""
    if shutil.which("bcl2fastq") is not None:
        return "bcl2fastq"
    elif shutil.which("bcl-convert") is not None:
        return "bcl-convert"  
    else:
        logging.error(
            "bcl2fastq or bcl-convert are not found in the system PATH. Please install bcl2fastq or bcl-convert and make sure it's in your PATH."
        )
    
    return None


def check_os():
    """Check if the operating system is Linux."""
    if sys.platform != "linux":
        logging.warning(
            "This script is designed to run on Linux. You may encounter issues on other operating systems."
        )
        return False
    return True

def create_lanes_tiles_S4():
    lanes_and_tiles = []
    for lane in range(1, 5):
        for side in range(1, 3):
            for column in range(1, 7):
                for row in range(1, 79):
                    lanes_and_tiles.append(f"{lane}_{side}{column}{row:02d}")

    return lanes_and_tiles

def _run_flowcell_map(args: argparse.Namespace):
    if not check_os():
        return

    demux_tool = check_bcl2fastq()
    if demux_tool is None:
        return
    
    args.demux_tool = demux_tool

    max_cores = multiprocessing.cpu_count()
    if args.parallel_processes > max_cores:
        logging.warning(
            f"Requested {args.parallel_processes} processes, but only {max_cores} cores are available. Limiting to {max_cores} processes."
        )
        args.parallel_processes = max_cores

    if not os.path.exists(args.tiles_out):
        os.makedirs(args.tiles_out)
        logging.warning(f"Created output directory: {args.tiles_out}")

    run_info_path = os.path.join(args.bcl_in, "RunInfo.xml")
    if not os.path.exists(run_info_path):
        lanes_and_tiles = parse_run_info_xml(run_info_path)
    else:
        logging.warning(f"RunInfo.xml not found at {run_info_path}. Generating lanes and tiles programmatically.")
        lanes_and_tiles = create_lanes_tiles_S4()
        

    with tempfile.TemporaryDirectory(dir=args.tiles_out) as temp_dir:
        args.bcl_out = os.path.join(temp_dir, "bcl_out")
        args.tilecoords_out = os.path.join(temp_dir, "tilecoords_out")
        args.dedup_out = os.path.join(temp_dir, "dedup_out")
        args.merge_out = os.path.join(temp_dir, "merge_out")
        args.distribute_out = os.path.join(temp_dir, "distribute_out")

        for dir_path in [args.bcl_out, args.tilecoords_out, args.dedup_out, args.merge_out, args.distribute_out]:
            os.makedirs(dir_path, exist_ok=True)
    
        lanes_and_tiles_path = os.path.join(args.bcl_in, "lanes_and_tiles.txt")
        with open(lanes_and_tiles_path, "w") as f:
            for tile in lanes_and_tiles:
                f.write(f"{tile}\n")

        # process tiles in parallel
        with ProcessPoolExecutor(max_workers=args.parallel_processes) as executor:
            futures = [executor.submit(process_tile, tile, args) for tile in lanes_and_tiles]
            for _ in tqdm(futures, total=len(futures), desc="Processing tiles"):
                _.result()

        if not os.path.exists(args.tilecoords_out) or not os.listdir(args.tilecoords_out):
            logging.error(f"Tile coordinates output directory {args.tilecoords_out} does not exist or is empty")
            return

        logging.info("Deduplicating individual tiles")
        deduplicate_tiles(args.tilecoords_out, args.dedup_out)

        if not os.path.exists(args.dedup_out) or not os.listdir(args.dedup_out):
            logging.error(f"Deduplicated tiles directory {args.dedup_out} does not exist or is empty")
            return

        logging.info("Merging and deduplicating all tiles")
        merged_file = os.path.join(args.merge_out, "merged_deduplicated.txt")
        merge_tiles(args.dedup_out, merged_file, args.parallel_processes)

        if not os.path.exists(merged_file):
            logging.error(f"Merged file {merged_file} does not exist")
            return

        logging.info("Writing tiles")
        distribute_per_file(merged_file, args.distribute_out)

        if not os.path.exists(args.distribute_out) or not os.listdir(args.distribute_out):
            logging.error(f"Distributed files directory {args.distribute_out} does not exist or is empty")
            return

        logging.info("Compressing tiles")
        compress_files(args.distribute_out, args.tiles_out)

    if not os.path.exists(args.tiles_out) or not os.listdir(args.tiles_out):
        logging.error(f"Final puck file directory {args.tiles_out} does not exist or is empty")
        return

    # TODO: create a file with statistics (how many tiles are created, how many barcodes, how many deduplicated)
    logging.info("Flowcell mapping completed successfully")


if __name__ == "__main__":
    from openst.cli import get_flowcell_map_parser

    args = get_flowcell_map_parser().parse_args()
    _run_flowcell_map(args)
