import logging

import numpy as np
import pandas as pd
from anndata import AnnData
from openst.utils.scanpy.pp import calculate_qc_metrics
from scipy.sparse import csc_matrix, csr_matrix, dok_matrix, vstack
from tqdm import tqdm


def calculate_adata_metrics(adata, dge_summary_path=None, n_reads=None):
    # calculate mitochondrial gene percentage
    adata.var["mt"] = (
        adata.var_names.str.startswith("Mt-")
        | adata.var_names.str.startswith("mt-")
        | adata.var_names.str.startswith("MT-")
    )

    calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    add_reads = False
    if dge_summary_path is not None:
        dge_summary = pd.read_csv(
            dge_summary_path,
            skiprows=7,
            sep="\t",
            index_col="cell_bc",
            names=["cell_bc", "n_reads", "n_umi", "n_genes"],
        )

        adata.obs = pd.merge(adata.obs, dge_summary[["n_reads"]], left_index=True, right_index=True)

        add_reads = True

    if n_reads is not None:
        adata.obs["n_reads"] = n_reads
        add_reads = True

    if add_reads:
        adata.obs["reads_per_counts"] = adata.obs['n_reads'] / adata.obs['total_counts']


def reassign_indices_adata(adata, new_ilocs, joined_coordinates, labels):
    original_ilocs = np.arange(new_ilocs.shape[0])
    sorted_ix = np.argsort(new_ilocs)
    new_ilocs = new_ilocs[sorted_ix]
    original_ilocs = original_ilocs[sorted_ix]

    joined_C = adata.X[original_ilocs]

    change_ix = np.where(new_ilocs[:-1] != new_ilocs[1:])[0] + 1

    ix_array = np.asarray(np.split(np.arange(new_ilocs.shape[0]), change_ix, axis=0), dtype="object")

    joined_C_sumed = vstack(
        [csr_matrix(joined_C[ix_array[n].astype(int), :].sum(0)) for n in tqdm(range(len(ix_array)), desc="Summarising expression per segmented cell")]
    )

    adata_out = AnnData(
        csc_matrix(joined_C_sumed),
        obs=pd.DataFrame(
            {"x_pos": joined_coordinates[:, 0], "y_pos": joined_coordinates[:, 1], "cell_ID_mask": labels}
        ),
        var=adata.var,
    )

    adata_out.obsm["spatial"] = joined_coordinates

    # rename index
    adata_out.obs.index.name = "cell_bc"

    def summarise_adata_obs_column(adata, column, summary_fun=sum):
        vals_to_join = adata.obs[column].to_numpy()[original_ilocs]
        vals_joined = np.array(
            [summary_fun(vals_to_join[ix_array[n].astype(int)]) for n in tqdm(range(len(ix_array)), desc=f"Summarising {column}")]
        )
        return vals_joined

    # summarise and attach n_reads, calculate metrics (incl. pcr)
    calculate_adata_metrics(
        adata_out,
        # provide the n_reads as a parameter
        n_reads=summarise_adata_obs_column(adata, "n_reads"),
    )

    adata_out.obs["n_joined"] = [len(x) for x in ix_array]

    mesh_bc_ilocs = np.arange(len(original_ilocs))[original_ilocs]

    joined_dict = {i: mesh_bc_ilocs[ix_array[i].astype(int)] for i, x in enumerate(ix_array)}

    indices_joined_spatial_units = dok_matrix((len(joined_dict), len(adata.obs_names)), dtype=np.int8)

    for obs_name_aggregate, obs_name_to_aggregate in joined_dict.items():
        indices_joined_spatial_units[obs_name_aggregate, obs_name_to_aggregate] = 1

    indices_joined_spatial_units = indices_joined_spatial_units.tocsr()
    adata_out.uns["spatial_units_obs_names"] = np.array(adata.obs_names)
    adata_out.uns["indices_joined_spatial_units"] = indices_joined_spatial_units

    from statistics import mean

    for column in ["exact_entropy", "theoretical_entropy", "exact_compression", "theoretical_compression"]:
        adata_out.obs[column] = summarise_adata_obs_column(adata, column, mean)

    return adata_out
