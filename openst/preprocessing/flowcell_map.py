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

import csv
import gzip
import resource
import shutil

from tqdm import tqdm

# for the last merge sort, in gigabytes
MEM_PER_CORE = 4
CHUNK_SIZE = 10_000_000
MAX_FILES = 50_000

log_lock = threading.Lock()

def calculate_offsets(file_path, num_processes):
    file_size = os.path.getsize(file_path)
    chunk_size = file_size // num_processes
    offsets = []

    with open(file_path, 'rb') as f:
        offset = 0
        for _ in range(num_processes):
            offsets.append(offset)
            f.seek(chunk_size, 1)
            f.readline()  # Move to the next full line
            offset = f.tell()

    return offsets

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

def process_fastq(fq: str, args: argparse.Namespace):
    os.makedirs(args.tilecoords_out, exist_ok=True)

    openst_cmd = [
        "openst",
        "barcode_preprocessing",
        "--fastq-in",
        fq,
        "--tilecoords-out",
        args.tilecoords_out,
        "--out-suffix",
        args.out_suffix,
        "--out-prefix",
        args.out_prefix,
        "--crop-seq",
        args.crop_seq
    ]
    if args.rev_comp:
        openst_cmd.append("--rev-comp")

    returncode, stdout, stderr = run_command(openst_cmd, f"openst barcode_preprocessing")
    log_output(f"openst barcode_preprocessing for a multi-tile fastq file", returncode, stdout, stderr)
    if returncode != 0:
        logging.warning(f"return code {returncode}, at command '{openst_cmd}'")

def process_deduplication_file(file: str, input_folder: str, output_folder: str):
    input_file = os.path.join(input_folder, file)
    output_file = os.path.join(output_folder, f".uncompressed.sorted.{file}")
    cmd = f"""
    LC_ALL=C zcat "{input_file}" | 
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

def sort_chunk(chunk_files, output_file, temp_dir):
    cmd = f"LC_ALL=C sort -m -u -k1,1 -t $'\\t' -T {temp_dir} {' '.join(chunk_files)} > {output_file}"
    subprocess.run(["bash", "-c", cmd], check=True)
    return output_file

def merge_tiles(input_folder: str, output_file: str, threads: int, chunk_size: int = 100):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    files = [os.path.join(input_folder, f) for f in os.listdir(input_folder) if f.startswith(".uncompressed.sorted.")]
    
    # Create a temporary directory for intermediate files
    with tempfile.TemporaryDirectory(dir=input_folder) as temp_dir:
        # Stage 1: Sort chunks of files in parallel
        chunk_outputs = []
        with ProcessPoolExecutor(max_workers=threads) as executor:
            futures = []
            for i in range(0, len(files), chunk_size):
                chunk = files[i:i+chunk_size]
                chunk_output = os.path.join(temp_dir, f"chunk_sorted_{i}.txt")
                futures.append(executor.submit(sort_chunk, chunk, chunk_output, temp_dir))
            
            for future in tqdm(futures, desc="Sorting chunks"):
                chunk_outputs.append(future.result())
        
        # Stage 2: Final merge of all chunk outputs
        final_merge_cmd = f"LC_ALL=C sort --parallel={threads} -m -u -k1,1 -t $'\\t' -T {temp_dir} {' '.join(chunk_outputs)} > {output_file}"
        
        try:
            subprocess.run(["bash", "-c", final_merge_cmd], check=True)
            logging.info(f"Final merge completed successfully. Output: {output_file}")
        except subprocess.CalledProcessError as e:
            logging.error(f"Error during final merge: {str(e)}")
            raise


def process_chunk(file_path, start_offset, end_offset, output_folder, chunk_id, num_processes):
    file_handles = {}
    csvw_handles = {}
    temp_folder = os.path.join(output_folder, f"temp_{chunk_id}")
    os.makedirs(temp_folder, exist_ok=True)

    flush_interval = 10_000_000
    lines_processed = 0

    with open(file_path, 'r') as f:
        f.seek(start_offset)
        current_position = start_offset
        
        # Skip the first line if not starting from the beginning (it might be partial)
        if start_offset > 0:
            f.readline()
            current_position = f.tell()

        pbar = None
        if start_offset == 0:
            total_size = end_offset - start_offset
            pbar = tqdm(total=total_size*num_processes, unit='B', unit_scale=True, desc=f'Processing in {num_processes} chunks')
        
        while current_position < end_offset:
            line = f.readline()
            if not line:  # End of file
                break
            
            new_position = f.tell()
            if pbar:
                pbar.update((new_position - current_position)*num_processes)
            current_position = new_position

            row = next(csv.reader([line], delimiter='\t'))
            
            id = row[-1]
            if id not in file_handles:
                file_handles[id] = open(os.path.join(temp_folder, f"{id}"), 'w', newline='')
                csvw_handles[id] = csv.writer(file_handles[id], delimiter='\t')
                csvw_handles[id].writerow(["cell_bc", "xcoord", "ycoord"])
            csvw_handles[id].writerow(row[:-1])

            lines_processed += 1
            if lines_processed % flush_interval == 0:
                for handle in file_handles.values():
                    handle.flush()
                    os.fsync(handle.fileno())

        if pbar:
            pbar.close()

    for handle in file_handles.values():
        handle.flush()
        os.fsync(handle.fileno())
        handle.close()

def merge_file(base_name, temp_folders, output_folder):
    final_file = os.path.join(output_folder, f"{base_name}.txt.gz")
    temp_files = [os.path.join(folder, f"{base_name}.txt.gz") 
                  for folder in temp_folders 
                  if os.path.exists(os.path.join(folder, f"{base_name}.txt.gz"))]
    
    with gzip.open(final_file, 'wt', newline='') as f_out:
        writer = csv.writer(f_out, delimiter='\t')
        writer.writerow(["cell_bc", "xcoord", "ycoord"])  # Write header once
        for temp_file in temp_files:
            with open(temp_file, 'r') as f_in:
                reader = csv.reader(f_in, delimiter='\t')
                next(reader)  # Skip header
                for row in reader:
                    writer.writerow(row)
            os.remove(temp_file)
    
    return base_name

def merge_intermediate_files(output_folder, num_processes):
    temp_folders = [os.path.join(output_folder, f"temp_{i}") for i in range(num_processes)]
    
    # Get unique base names across all temp folders
    base_names = set()
    for folder in temp_folders:
        if not os.path.exists(folder):
            logging.error(f"The temporary folder '{folder}' does not exist")
            break
        base_names.update(name.rsplit('.', 2)[0] for name in os.listdir(folder))
    
    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        futures = [executor.submit(merge_file, base_name, temp_folders, output_folder) 
                   for base_name in base_names]
        
        for future in tqdm(futures, total=len(futures), desc="Merging files"):
            future.result()
    
    # Clean up temp folders
    for folder in temp_folders:
        shutil.rmtree(folder)

def distribute_per_file(input_file: str, output_folder: str, num_processes: int = None):
    os.makedirs(output_folder, exist_ok=True)
    
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    resource.setrlimit(resource.RLIMIT_NOFILE, (min(MAX_FILES, hard), hard))
    
    if num_processes is None:
        num_processes = os.cpu_count()
    
    try:
        offsets = calculate_offsets(input_file, num_processes)
        offsets.append(os.path.getsize(input_file))  # Add end of file offset
        
        with ProcessPoolExecutor(max_workers=num_processes) as executor:
            futures = []
            for i in range(num_processes):
                start_offset = offsets[i]
                end_offset = offsets[i+1]
                futures.append(executor.submit(process_chunk, input_file, start_offset, end_offset, output_folder, i, num_processes))
            
            for future in tqdm(futures, total=num_processes, desc="Processing chunks"):
                future.result()
        
        logging.info("All chunks processed. Merging intermediate files...")
        merge_intermediate_files(output_folder, num_processes)
    
    except Exception as e:
        logging.error(f"Error during file distribution: {str(e)}")
        raise
    
    logging.info(f"File distribution completed. Output in {output_folder}")


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

    if len(args.bcl_in) == 0 and len(args.fastq_in) == 0:
        logging.error("You must specify either --bcl-in <folder_to_bcl> or --fastq-in <fastq_files>")
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

    if os.path.exists(os.path.join(args.tiles_out, "done")):
        logging.warning("openst flowcell_map was executed successfully on this folder")
        return

    if not os.path.exists(args.tiles_out):
        os.makedirs(args.tiles_out)
        logging.warning(f"Created output directory: {args.tiles_out}")

    run_info_path = os.path.join(args.bcl_in, "RunInfo.xml")
    if os.path.exists(run_info_path):
        lanes_and_tiles = parse_run_info_xml(run_info_path)
    else:
        logging.warning(f"RunInfo.xml not found at {run_info_path}. Generating lanes and tiles programmatically.")
        lanes_and_tiles = create_lanes_tiles_S4()

    if args.tmp_dir == "" or args.tmp_dir is None:
        logging.debug("Will write temporary files to the same output directory.")
        args.tmp_dir = args.tiles_out

        
    args.bcl_out = os.path.join(args.tmp_dir, "bcl_out")
    args.tilecoords_out = os.path.join(args.tmp_dir, "tilecoords_out")
    args.dedup_out = os.path.join(args.tmp_dir, "dedup_out")
    args.merge_out = os.path.join(args.tmp_dir, "merge_out")
    args.distribute_out = os.path.join(args.tmp_dir, "distribute_out")
    merged_file = os.path.join(args.merge_out, "merged_deduplicated.txt")

    for dir_path in [args.bcl_out, args.tilecoords_out, args.dedup_out, args.merge_out, args.distribute_out]:
        os.makedirs(dir_path, exist_ok=True)

    lanes_and_tiles_path = os.path.join(args.tiles_out, "lanes_and_tiles.txt")
    with open(lanes_and_tiles_path, "w") as f:
        for tile in lanes_and_tiles:
            f.write(f"{tile}\n")

    if len(os.listdir(args.dedup_out)) == 0:
        # process tiles in parallel
        if len(args.bcl_in) != 0:
            with ProcessPoolExecutor(max_workers=args.parallel_processes) as executor:
                futures = [executor.submit(process_tile, tile, args) for tile in lanes_and_tiles]
                for _ in tqdm(futures, total=len(futures), desc="Processing tiles"):
                    _.result()
        elif len(args.fastq_in) != 0:
            with ProcessPoolExecutor(max_workers=args.parallel_processes) as executor:
                futures = [executor.submit(process_fastq, fq, args) for fq in args.fastq_in]
                for _ in tqdm(futures, total=len(futures), desc="Processing FASTQ files"):
                    _.result()

    if not os.path.exists(args.tilecoords_out) or not os.listdir(args.tilecoords_out):
        logging.error(f"Tile coordinates output directory {args.tilecoords_out} does not exist or is empty")
        return

    logging.info("Deduplicating individual tiles")
    if not os.path.exists(merged_file):
        deduplicate_tiles(args.tilecoords_out, args.dedup_out)

    if not os.path.exists(args.dedup_out) or not os.listdir(args.dedup_out):
        logging.error(f"Deduplicated tiles directory {args.dedup_out} does not exist or is empty")
        return

    logging.info("Merging and deduplicating all tiles")
    if not os.path.exists(os.path.join(args.merge_out, "done")):
        merge_tiles(args.dedup_out, merged_file, args.parallel_processes)

    open(os.path.join(args.merge_out, "done"), 'a').close()

    if not os.path.exists(merged_file):
        logging.error(f"Merged file {merged_file} does not exist")
        return
    
    logging.info("Writing compressed tiles")
    if not os.path.exists(os.path.join(args.tiles_out, "done")):
        distribute_per_file(merged_file, args.tiles_out, args.parallel_processes)

    open(os.path.join(args.tiles_out, "done"), 'a').close()

    if not os.path.exists(args.tiles_out) or not os.listdir(args.tiles_out):
        logging.error(f"Final puck file directory {args.tiles_out} does not exist or is empty")
        return
    
    logging.info("Cleaning temporary directories and files")
    for d in [args.bcl_out, args.tilecoords_out, args.dedup_out, args.merge_out, args.distribute_out]:
        if os.path.exists(d):
            os.remove(d)

    # TODO: create a file with statistics (how many tiles are created, how many barcodes, how many deduplicated)
    logging.info("Flowcell mapping completed successfully")


if __name__ == "__main__":
    from openst.cli import get_flowcell_map_parser

    args = get_flowcell_map_parser().parse_args()
    _run_flowcell_map(args)
