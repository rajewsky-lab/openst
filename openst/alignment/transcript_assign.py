import numpy as np
import pandas as pd
import scanpy as sc

from PIL import Image
from skimage import measure
from tqdm import tqdm

from skimage.segmentation import expand_labels
from scipy.sparse import csr_matrix, csc_matrix, vstack, dok_matrix

Image.MAX_IMAGE_PIXELS = 933120000

def setup_parser(parser):
    parser.add_argument(
        "--adata",
        type=str,
        help="path to previously aligned spatial.h5ad AnnData file",
        required=True,
    )

    parser.add_argument(
        "--mask",
        type=str,
        help="path to image mask; must be in same coordinates as the obsm['spatial'] in the AnnData file",
        required=True,
    )

    parser.add_argument(
        "--output",
        type=str,
        help="path and filename for output file that will be generated",
        required=True,
    )

    parser.add_argument(
        "--dilate-px",
        type=int,
        help="specify how many pixels the outlines of the segmentation mask will be extended",
        required=False,
        default=0
    )

    return parser

def calculate_adata_metrics(adata, dge_summary_path=None, n_reads=None):
    import scanpy as sc
    import pandas as pd

    # calculate mitochondrial gene percentage
    adata.var["mt"] = (
        adata.var_names.str.startswith("Mt-")
        | adata.var_names.str.startswith("mt-")
        | adata.var_names.str.startswith("MT-")
    )

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    add_reads = False
    if dge_summary_path is not None:
        dge_summary = pd.read_csv(
            dge_summary_path,
            skiprows=7,
            sep="\t",
            index_col="cell_bc",
            names=["cell_bc", "n_reads", "n_umi", "n_genes"],
        )

        adata.obs = pd.merge(
            adata.obs, dge_summary[["n_reads"]], left_index=True, right_index=True
        )

        add_reads = True

    if n_reads is not None:
        adata.obs["n_reads"] = n_reads
        add_reads = True

    if add_reads:
        adata.obs["reads_per_counts"] = adata.obs.n_reads / adata.obs.total_counts

def reassign_indices_adata(adata, new_ilocs, joined_coordinates, mask_image, labels):
    import pandas as pd
    import scanpy as sc
    import numpy as np
    import anndata
    
    original_ilocs = np.arange(new_ilocs.shape[0])
    sorted_ix = np.argsort(new_ilocs)
    new_ilocs = new_ilocs[sorted_ix]
    original_ilocs = original_ilocs[sorted_ix]

    joined_C = adata.X[original_ilocs]

    change_ix = np.where(new_ilocs[:-1] != new_ilocs[1:])[0] + 1

    ix_array = np.asarray(np.split(np.arange(new_ilocs.shape[0]), change_ix, axis=0), dtype='object')

    joined_C_sumed = vstack([csr_matrix(joined_C[ix_array[n].astype(int), :].sum(0)) for n in tqdm(range(len(ix_array)))])

    adata_out = anndata.AnnData(csc_matrix(joined_C_sumed), 
        obs = pd.DataFrame({'x_pos': joined_coordinates[:, 0],
                            'y_pos': joined_coordinates[:, 1],
                            'cell_ID_mask': labels}),
        var = adata.var)

    adata_out.obsm['spatial'] = joined_coordinates
    adata_out.uns['spatial'] = {'mask': mask_image}

    # rename index
    adata_out.obs.index.name = 'cell_bc'

    def summarise_adata_obs_column(adata, column, summary_fun=sum):
        vals_to_join = adata.obs[column].to_numpy()[original_ilocs]
        vals_joined = np.array(
            [summary_fun(vals_to_join[ix_array[n].astype(int)])
                for n in tqdm(range(len(ix_array)))])
        return vals_joined

    # summarise and attach n_reads, calculate metrics (incl. pcr)
    calculate_adata_metrics(adata_out,
        # provide the n_reads as a parameter
        n_reads = summarise_adata_obs_column(adata, 'n_reads'))

    adata_out.obs['n_joined'] = [len(x) for x in ix_array]

    mesh_bc_ilocs = np.arange(len(original_ilocs))[original_ilocs]

    joined_dict = {i: mesh_bc_ilocs[x] for i, x in enumerate(ix_array)}

    indices_joined_spatial_units = dok_matrix(
        (len(joined_dict), len(adata.obs_names)), dtype=np.int8
    )

    for obs_name_aggregate, obs_name_to_aggregate in joined_dict.items():
        indices_joined_spatial_units[obs_name_aggregate, obs_name_to_aggregate] = 1

    indices_joined_spatial_units = indices_joined_spatial_units.tocsr()
    adata_out.uns["spatial_units_obs_names"] = np.array(adata.obs_names)
    adata_out.uns["indices_joined_spatial_units"] = indices_joined_spatial_units

    from statistics import mean

    for column in ['exact_entropy', 'theoretical_entropy', 'exact_compression',\
        'theoretical_compression']:
        adata_out.obs[column] = summarise_adata_obs_column(adata, column, mean)

    return adata_out

def transfer_segmentation(adata_transformed_coords, label_image, props_filter):
    joined_coordinates = np.array([props_filter['centroid-0'], props_filter['centroid-1']]).T
    joined_coordinates = np.vstack([np.array([0, 0]), joined_coordinates])

    cell_ID_merged = np.array(props_filter['label'])
    cell_ID_merged = np.hstack([np.array([0]), cell_ID_merged])

    adata_by_cell = reassign_indices_adata(adata_transformed_coords, np.array(adata_transformed_coords.obs["cell_ID"]), joined_coordinates, label_image, cell_ID_merged)

    spatial_units_obs_names_dict = {}

    for sn in tqdm(adata_by_cell.uns['spatial_units_obs_names']):
        bc, tile = sn.split(":")
        if tile in spatial_units_obs_names_dict.keys():
            spatial_units_obs_names_dict[tile] += [bc]
        else:
            spatial_units_obs_names_dict[tile] = [bc]

    for sn in tqdm(adata_by_cell.uns['spatial_units_obs_names']):
        bc, tile = sn.split(":")
        if tile in spatial_units_obs_names_dict.keys():
            spatial_units_obs_names_dict[tile] += [bc]
        else:
            spatial_units_obs_names_dict[tile] = [bc]

    return adata_by_cell

def prepare_mask_and_adata(mask, adata, dilate_px : int = 0 ):
    if dilate_px > 0:
        _mask = expand_labels(mask, distance=dilate_px)
    else:
        _mask = mask

    label_image = measure.label(_mask)
    adata = adata[(adata.obsm['spatial'][:, 0] <= label_image.shape[0]) &
                                                        (adata.obsm['spatial'][:, 1] <= label_image.shape[1])]

    labels = label_image[adata.obsm['spatial'][:, 0].astype(int),
                        adata.obsm['spatial'][:, 1].astype(int)]

    adata.obs["cell_ID"] = labels

    props = measure.regionprops_table(label_image,
                            properties=['label', 'centroid', 'area', 'axis_major_length', 'axis_minor_length'])
    props = pd.DataFrame(props)

    props_filter = props[props.label.isin(np.unique(adata.obs["cell_ID"]))]
    return _mask, adata, props_filter

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description="transfer a segmentation mask to the individual spots from a previously aligned AnnData file",
    )

    parser = setup_parser(parser)
    args = parser.parse_args()

    print("Loading aligned AnnData")
    adata_transformed_coords = sc.read_h5ad(args.adata)
    non_barcodes = pd.Series(adata_transformed_coords.obs_names).apply(lambda x: x.split(":")[0].isnumeric()).values
    adata_transformed_coords = adata_transformed_coords[~non_barcodes].copy()

    print("Loading image mask")
    stitched_image_he_mask = np.array(Image.open(args.mask))
    
    print("Preparing segmentation mask and fixing AnnData spatial coordinates")
    label_image, adata_transformed_coords, props_filter = prepare_mask_and_adata(stitched_image_he_mask, adata_transformed_coords, args.dilate_px)

    print("Transferring spatial data to segmentation mask")
    adata_by_cell = transfer_segmentation(adata_transformed_coords, label_image, props_filter)

    print("Writing output")
    adata_by_cell.write_h5ad(args.output)