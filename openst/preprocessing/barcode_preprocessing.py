import gzip
import os
import time

import pandas as pd

import logging
import dnaio

tab = str.maketrans("ACTG", "TGAC")

def reverse_complement_table(seq):
    return seq.translate(tab)[::-1]


def get_tile_number_and_coordinates(read_name):
    read_name = read_name.split(" ")[0].split(":")[-4:]
    return read_name


def process_multiple_tiles(
    in_fastq: str, out_path: str, out_prefix: str, out_suffix: str, sequence_preprocessor: callable = None
):
    current_tile_number = None
    n_written = 0
    tile_file = None

    def open_output_file(filename: str):
        """Helper function to open file with appropriate compression"""
        if filename.endswith('.gz'):
            return gzip.open(filename, 'wt')
        return open(filename, 'w')

    try:
        with dnaio.open(in_fastq) as f:
            for line in f:
                lane, tile_number, xcoord, ycoord = get_tile_number_and_coordinates(line.name)
                if current_tile_number != tile_number:

                    if current_tile_number is not None:
                        logging.info(f"Written {n_written} spatial barcodes to {lane}_{current_tile_number}")
                        tile_file.close()

                    current_tile_number = tile_number

                    tile_fname = os.path.join(
                                    out_path,
                                    out_prefix + lane + "_" + current_tile_number + out_suffix,
                                )

                    tile_file = open_output_file(tile_fname)
                    tile_file.write("\t".join(["cell_bc", "xcoord", "ycoord"])+"\n")
                    n_written = 0

                sequence = sequence_preprocessor(line.sequence) if sequence_preprocessor is not None else line.sequence
                tile_file.write("\t".join([sequence, xcoord, ycoord])+"\n")
                n_written += 1
    finally:
        if tile_file is not None:
            if current_tile_number is not None:
                logging.info(f"Written {n_written} spatial barcodes to {lane}_{current_tile_number}")
            tile_file.close()

def process_single_tile(in_fastq: str, sequence_preprocessor: callable = None) -> pd.DataFrame:
    all_tile_numbers = set()
    sequences, xcoords, ycoords = [[], [], []]
    idx = 0
    with gzip.open(in_fastq, "rt") as f:
        for line in f:
            if idx % 4 == 0:
                _, tile_number, xcoord, ycoord = get_tile_number_and_coordinates(line.strip())
                all_tile_numbers.add(tile_number)
                if len(all_tile_numbers) > 1:
                    raise ValueError("You specified a single tile, and the file contains more than one...")

            if idx % 4 == 1:
                sequence = sequence_preprocessor(line) if sequence_preprocessor is not None else line
                sequences.append(sequence)
                xcoords.append(xcoord)
                ycoords.append(ycoord)

            idx += 1

    df = pd.DataFrame({"cell_bc": sequences, "x_pos": xcoords, "y_pos": ycoords})
    return df


def _run_barcode_preprocessing(args):
    crop_seq_slice = slice(
        *[{True: lambda n: None, False: int}[x == ""](x) for x in (args.crop_seq.split(":") + ["", "", ""])[:3]]
    )

    def sequence_preprocessor(sequence):
        sequence = sequence[crop_seq_slice].strip()
        if args.rev_comp:
            sequence = reverse_complement_table(sequence)
        return sequence

    start_time = time.time()

    if not args.single_tile and not args.unsorted:
        process_multiple_tiles(args.fastq_in, args.tilecoords_out, args.out_prefix, args.out_suffix, sequence_preprocessor)
    elif not args.single_tile and args.unsorted:
        raise NotImplementedError("We don't support multi-tile files with unsorted tiles yet!")
    else:
        df = process_single_tile(args.fastq_in, sequence_preprocessor)
        logging.info(f"Writing {len(df):,} barcodes to {os.path.join(args.tilecoords_out, args.out_prefix + args.out_suffix)}")
        df.to_csv(
            os.path.join(args.tilecoords_out, args.out_prefix + args.out_suffix),
            index=False,
            sep="\t",
        )

    logging.info(f"Finished in {round(time.time()-start_time, 2)} sec")


if __name__ == "__main__":
    from openst.cli import get_barcode_preprocessing_parser
    args = get_barcode_preprocessing_parser().parse_args()
    _run_barcode_preprocessing(args)
