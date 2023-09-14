import argparse
import gzip
import os
import time

import pandas as pd


def get_barcode_preprocessing_parser():
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        add_help=False,
        description="convert fastq files into spatial barcode files (sequence and coordinates)",
    )
    parser.add_argument("--in-fastq", type=str, required=True, help="path to the fastq file")
    parser.add_argument(
        "--out-path",
        type=str,
        required=True,
        help="folder where the output files will be generated",
    )
    parser.add_argument(
        "--out-suffix",
        type=str,
        required=True,
        help="where to write the output file. it works as the suffix when multiple tiles are generated",
    )
    parser.add_argument(
        "--out-prefix",
        type=str,
        default="",
        help="where to write the output file. it works as the prefix when multiple tiles are generated",
    )
    parser.add_argument(
        "--crop-seq",
        type=str,
        default=":",
        help="crop the input sequence, should be a valid python slice",
    )
    parser.add_argument(
        "--rev-comp",
        action="store_true",
        help="applies reverse complementary after sequence cropping",
    )
    parser.add_argument(
        "--single-tile",
        action="store_true",
        help="it is guarranteed that the input .fastq(.gz) file contains only a tile. Throw an error otherwise",
    )
    parser.add_argument(
        "--unsorted",
        action="store_true",
        help="supports that the file is unsorted respect to tiles. might be slower",
    )

    return parser


def setup_barcode_preprocessing_parser(parent_parser):
    """setup_barcode_preprocessing_parser"""
    parser = parent_parser.add_parser(
        "barcode_preprocessing",
        help="convert fastq files into spatial barcode files (sequence and coordinates)",
        parents=[get_barcode_preprocessing_parser()],
    )
    parser.set_defaults(func=_run_barcode_preprocessing)

    return parser


tab = str.maketrans("ACTG", "TGAC")


def reverse_complement_table(seq):
    return seq.translate(tab)[::-1]


def get_tile_number_and_coordinates(read_name):
    read_name = read_name.split(" ")[0].split(":")[-3:]
    return read_name


def process_multiple_tiles(
    in_fastq: str, out_path: str, out_prefix: str, out_suffix: str, sequence_preprocessor: callable = None
):
    current_tile_number = None
    sequences, xcoords, ycoords = [[], [], []]
    idx = 0
    with gzip.open(in_fastq, "rt") as f:
        for line in f:
            if idx % 4 == 0:
                tile_number, xcoord, ycoord = get_tile_number_and_coordinates(line.strip())
                if current_tile_number is None:
                    current_tile_number = tile_number

                if tile_number != current_tile_number:
                    print(f"Writing {len(sequences):,} barcodes of file {current_tile_number} to disk")
                    df = pd.DataFrame({"cell_bc": sequences, "xcoord": xcoords, "ycoord": ycoords})
                    df.to_csv(
                        os.path.join(
                            out_path,
                            out_prefix + current_tile_number + out_suffix,
                        ),
                        index=False,
                        sep="\t",
                    )
                    current_tile_number = tile_number
                    sequences, xcoords, ycoords = [[], [], []]
            elif idx % 4 == 1:
                sequence = sequence_preprocessor(line) if sequence_preprocessor is not None else line
                sequences.append(sequence)
                xcoords.append(xcoord)
                ycoords.append(ycoord)

            idx += 1


def process_single_tile(in_fastq: str, sequence_preprocessor: callable = None) -> pd.DataFrame:
    all_tile_numbers = set()
    sequences, xcoords, ycoords = [[], [], []]
    idx = 0
    with gzip.open(in_fastq, "rt") as f:
        for line in f:
            if idx % 4 == 0:
                tile_number, xcoord, ycoord = get_tile_number_and_coordinates(line.strip())
                all_tile_numbers.add(tile_number)
                if len(all_tile_numbers) > 1:
                    raise ValueError("You specified a single tile, and the file contains more than one...")

            if idx % 4 == 1:
                sequence = sequence_preprocessor(line) if sequence_preprocessor is not None else line
                sequences.append(sequence)
                xcoords.append(xcoord)
                ycoords.append(ycoord)

            idx += 1

    df = pd.DataFrame({"cell_bc": sequences, "xcoord": xcoords, "ycoord": ycoords})
    return df


def _run_barcode_preprocessing(args):
    crop_seq_slice = slice(
        *[{True: lambda n: None, False: int}[x == ""](x) for x in (args.crop_seq.split(":") + ["", "", ""])[:3]]
    )

    def sequence_preprocessor(sequence):
        sequence = sequence[crop_seq_slice].strip()
        if args.rev_comp:
            sequence = reverse_complement_table(sequence)

    start_time = time.time()

    if not args.single_tile and not args.unsorted:
        process_multiple_tiles(args.in_fastq, args.out_path, args.out_prefix, args.out_suffix, sequence_preprocessor)
    elif not args.single_tile and args.unsorted:
        raise NotImplementedError("We don't support multi-tile files with unsorted tiles yet!")
    else:
        df = process_single_tile(args.in_fastq, sequence_preprocessor)
        print(f"Writing {len(df):,} barcodes to {os.path.join(args.out_path, args.out_prefix + args.out_suffix)}")
        df.to_csv(
            os.path.join(args.out_path, args.out_prefix + args.out_suffix),
            index=False,
            sep="\t",
        )

    print(f"Finished in {round(time.time()-start_time, 2)} sec")


if __name__ == "__main__":
    args = get_barcode_preprocessing_parser().parse_args()
    _run_barcode_preprocessing(args)
