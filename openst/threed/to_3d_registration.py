import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd
from anndata import AnnData, read_h5ad

from openst.utils.file import check_directory_exists, check_file_exists


def convert_adata_to_crosstab(adata: AnnData, genes: list = None):
    if genes is not None and len(adata.var_names[adata.var_names.isin(genes)]) == 0:
        raise ValueError(f"The genes {genes} could not be found in the specified adata")

    if genes is not None:
        adata = adata[:, adata.var_names.isin(genes)]

    adata = adata[adata.obs_names[np.unique(adata.X.nonzero()[0])]].copy()

    df_locations = pd.DataFrame(
        {"barcodes": adata.obs_names, "x_pos": adata.obsm["spatial"][:, 0], "y_pos": adata.obsm["spatial"][:, 1]}
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

    # Check input and output data
    check_file_exists(args.in_adata)

    if not check_directory_exists(args.output_dir):
        raise FileNotFoundError("Parent directory for --output-mask does not exist")

    logging.info(f"Loading file from {args.in_adata}")
    adata = read_h5ad(args.in_adata)

    logging.info("Preprocessing adata")
    adata.obsm["spatial"] -= adata.obsm["spatial"].min(axis=0)
    if args.rescale != 1:
        adata.obsm["spatial"] = adata.obsm["spatial"] * args.rescale

    if args.filter_umi_min >= 0:
        # TODO: check if this dimension or other
        adata = adata[np.array(adata.X.sum(axis=0)).flatten() > args.filter_umi_min]
    if args.filter_umi_max != -1 and args.filter_umi_max > args.filter_umi_min:
        adata = adata[np.array(adata.X.sum(axis=0)).flatten() < args.filter_umi_max]

    logging.info("Converting adata to crosstab")
    df_locations, df_genes = convert_adata_to_crosstab(adata, args.genes)

    logging.info(f"Saving locations and genes into csv files at {args.output_dir}")
    in_adata_stem = Path(args.in_adata).stem
    df_locations.to_csv(os.path.join(args.output_dir, f"{in_adata_stem}.locations.csv"), index=False)
    df_genes.to_csv(os.path.join(args.output_dir, f"{in_adata_stem}.genes.csv"))


if __name__ == "__main__":
    from openst.cli import get_to_3d_registration_parser
    args = get_to_3d_registration_parser().parse_args()
    _run_to_3d_registration(args)
