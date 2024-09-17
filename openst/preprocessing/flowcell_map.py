import argparse
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor
from typing import List
import logging
from tqdm import tqdm
import xml.etree.ElementTree as ET
import multiprocessing

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

    # Run bcl2fastq
    bcl2fastq_cmd = [
        "bcl2fastq",
        "-R", args.bcl_in,
        "--no-lane-splitting",
        "-o", bcl_out_dir,
        f"--tiles s_{tile}"
    ]
    subprocess.run(bcl2fastq_cmd, check=True)

    # Run openst barcode_preprocessing
    openst_cmd = [
        "openst", "barcode_preprocessing",
        "--fastq-in", fastq_file,
        "--tilecoords-out", args.tilecoords_out,
        "--out-suffix", args.out_suffix,
        "--out-prefix", f"{args.out_prefix}{tile}",
        "--crop-seq", args.crop_seq,
        "--single-tile"
    ]
    if args.rev_comp:
        openst_cmd.append("--rev-comp")
    
    subprocess.run(openst_cmd, check=True)

def deduplicate_tiles(input_folder: str, output_folder: str):
    os.makedirs(output_folder, exist_ok=True)
    files = [f for f in os.listdir(input_folder) if f.endswith('.txt.gz')]

    def process_file(file: str):
        output_file = os.path.join(output_folder, f".uncompressed.sorted.{file}")
        subprocess.run(f"zcat {os.path.join(input_folder, file)} | tail -n +2 | sort -u -k1,1 -t$'\\t' > {output_file}", shell=True, check=True)

    with ProcessPoolExecutor() as executor:
        list(tqdm(executor.map(process_file, files), total=len(files), desc="Deduplicating tiles"))

def merge_tiles(input_folder: str, output_file: str):
    files = [os.path.join(input_folder, f) for f in os.listdir(input_folder) if f.startswith('.uncompressed.sorted.')]
    cmd = f"sort --batch-size=128 -m -u -k1,1 -t$'\\t' {' '.join(files)} > {output_file}"
    subprocess.run(cmd, shell=True, check=True)

def distribute_per_file(input_file: str, output_folder: str):
    os.makedirs(output_folder, exist_ok=True)
    cmd = f"awk -F'\\t' '{{print $1\"\\t\"$2\"\\t\"$3 > \"{output_folder}/\"$4}}' {input_file}"
    subprocess.run(cmd, shell=True, check=True)

def compress_files(input_folder: str, output_folder: str):
    os.makedirs(output_folder, exist_ok=True)
    files = [f for f in os.listdir(input_folder) if not f.startswith('.')]

    def process_file(file: str):
        input_file = os.path.join(input_folder, file)
        output_file = os.path.join(output_folder, f"{file}.gz")
        subprocess.run(f"gzip -c {input_file} > {output_file}", shell=True, check=True)

    with ProcessPoolExecutor() as executor:
        list(tqdm(executor.map(process_file, files), total=len(files), desc="Compressing files"))

def _run_flowcell_map(args: argparse.Namespace):
    max_cores = multiprocessing.cpu_count()
    if args.parallel_processes > max_cores:
        logging.warning(f"Requested {args.parallel_processes} processes, but only {max_cores} cores are available. Limiting to {max_cores} processes.")
        args.parallel_processes = max_cores

    run_info_path = os.path.join(args.bcl_in, "RunInfo.xml")
    if not os.path.exists(run_info_path):
        logging.error(f"RunInfo.xml not found at {run_info_path} - did you specify a proper Illumina basecalls folder?")
        return
    
    lanes_and_tiles = parse_run_info_xml(run_info_path)

    lanes_and_tiles_path = os.path.join(args.bcl_in, "lanes_and_tiles.txt")
    with open(lanes_and_tiles_path, 'w') as f:
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

    deduplicate_tiles(args.tilecoords_out, args.dedup_out)

    if not os.path.exists(args.dedup_out) or not os.listdir(args.dedup_out):
        logging.error(f"Deduplicated tiles directory {args.dedup_out} does not exist or is empty")
        return

    merged_file = os.path.join(args.merge_out, "merged_deduplicated.txt")
    merge_tiles(args.dedup_out, merged_file)

    if not os.path.exists(merged_file):
        logging.error(f"Merged file {merged_file} does not exist")
        return

    distribute_per_file(merged_file, args.distribute_out)

    if not os.path.exists(args.distribute_out) or not os.listdir(args.distribute_out):
        logging.error(f"Distributed files directory {args.distribute_out} does not exist or is empty")
        return

    compress_files(args.distribute_out, args.compress_out)

    if not os.path.exists(args.compress_out) or not os.listdir(args.compress_out):
        logging.error(f"Compressed files directory {args.compress_out} does not exist or is empty")
        return

    logging.info("Flowcell mapping completed successfully")

if __name__ == "__main__":
    from openst.cli import get_flowcell_map_parser
    args = get_flowcell_map_parser().parse_args()
    _run_flowcell_map(args)