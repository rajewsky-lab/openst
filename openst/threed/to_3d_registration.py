import argparse
import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from openst.utils.file import check_directory_exists, check_file_exists


def get_to_3d_registration_parser():
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="convert openst data for STIM (serial-section 3D registration); one file at a time",
        allow_abbrev=False,
        add_help=False,
    )
    parser.add_argument(
        "--in-adata",
        type=str,
        required=True,
        help="Path to the input adata file. Expects a file with raw counts.",
    )
    parser.add_argument(
        "--output-dir",
        type=int,
        default=-1,
        help="Path to the directory where the .genes and .locations csv files will be stored.",
    )
    parser.add_argument(
        "--filter-umi-max",
        type=int,
        default=-1,
        help="Preserve cells with at most < --filter-umi-max UMIs",
    )
    parser.add_argument(
        "--lognorm",
        action="store_true",
        help="When specified, applies scanpy's log normalization to the raw counts",
    )
    parser.add_argument(
        "--rescale",
        type=float,
        default=1,
        help="Rescales (multiples) the coordinates by --rescale units",
    )
    parser.add_argument(
        "--genes",
        type=str,
        nargs="+",
        deafult=None,
        help="Exports the expression level for the specified genes",
    )
    return parser


def setup_to_3d_registration_parser(parent_parser):
    """setup_to_3d_registration_parser"""
    parser = parent_parser.add_parser(
        "to_3d_registration",
        help="convert h5ad files to native STIM-supported format",
        parents=[get_to_3d_registration_parser()],
    )
    parser.set_defaults(func=_run_to_3d_registration)

    return parser


def convert_adata_to_crosstab(adata: AnnData, genes: list = None):
    if genes is not None and len(adata.var_names[adata.var_names.isin(genes)]) == 0:
        raise ValueError(f"The genes {genes} could not be found in the specified adata")

    if genes is not None:
        adata = adata[:, adata.var_names.isin(genes)]

    adata = adata[adata.obs_names[np.unique(adata.X.nonzero()[0])]].copy()

    df_locations = pd.DataFrame(
        {"barcodes": adata.obs_names, "xcoord": adata.obsm["spatial"][:, 0], "ycoord": adata.obsm["spatial"][:, 1]}
    )
    df_genes = pd.crosstab(
        np.repeat(adata.var_names[adata.X.nonzero()[1]], adata.X.data),
        np.repeat(adata.obs_names[adata.X.nonzero()[0]], adata.X.data),
    )

    return df_locations, df_genes


def _run_to_3d_registration(args):
    """
    _run_to_3d_registration
    """
    logging.info("openst to_3d_registrationation; running with parameters:")
    logging.info(args.__dict__)

    # Check input and output data
    check_file_exists(args.in_adata)

    if not check_directory_exists(args.output_dir):
        raise FileNotFoundError("Parent directory for --output-mask does not exist")

    logging.info(f"Loading file from {args.in_adata}")
    adata = sc.read_h5ad(args.in_adata)

    logging.info("Preprocessing adata")
    if args.filter_umi_min >= 0:
        sc.pp.filter_cells(adata, min_counts=args.filter_umi_min)
    if args.filter_umi_max != -1 and args.filter_umi_max > args.filter_umi_min:
        sc.pp.filter_genes(adata, max_counts=args.filter_umi_max)
    if args.lognorm:
        sc.pp.normalize_total(adata, inplace=True)
        sc.pp.log1p(adata)
    if args.rescale != 1:
        adata.obsm["spatial"] = adata.obsm["spatial"] * args.rescale

    logging.info("Converting adata to crosstab")
    df_locations, df_genes = convert_adata_to_crosstab(adata, args.genes)

    logging.info(f"Saving locations and genes into csv files at {args.output_dir}")
    in_adata_stem = Path(args.in_adata).stem
    df_locations.to_csv(os.path.join(args.output_dir, f"{in_adata_stem}.locations.csv"), index=False)
    df_genes.to_csv(os.path.join(args.output_dir, f"{in_adata_stem}.genes.csv"))


if __name__ == "__main__":
    args = get_to_3d_registration_parser().parse_args()
    _run_to_3d_registration(args)
