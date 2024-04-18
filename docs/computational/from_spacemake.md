# Running `openst` on `spacemake` folders
If you preprocess the sequencing data with `spacemake` (as we recommend), it is possible (since `openst>=0.0.7`)
to automatically populate data input/output arguments for the relevant commands (e.g., `segment`, `spatial_stitch`...).
This is done with `openst from_spacemake`, which the universal structure:

```sh
openst from_spacemake \
     --project-id <project_id> \
     --sample-id <sample_id> \
     --run-mode <one_of_run_modes> # optional, in case you specified >1 run_modes in the sample
     <subcommand> # it is one of image_stitch, spatial_stitch, pairwise_aligner
               # apply_transform, manual_pairwise_aligner, segment, segment_merge,
               # transcript_assign, merge_modalities, pseudoimage, preview
     <arguments_subcommand>
```

You can get help for the subcommands by running `openst <subcommand> --help`, as `openst from_spacemake <subcommand> --help` will only show the help for the `from_spacemake` command, not for `<subcommand>`. Most arguments displayed as `required` for the `<subcommand>` are overrided by `from_spacemake`.

## Getting started
To start, you need to `cd` into the root spacemake folder (i.e., containing the `project_df.csv` and `config.yaml` files).
You can check by `ls` inside that specific directory, you should get something like:

```sh
(spacemake) user@computer:~$ cd /openst_adult_demo/spacemake
(spacemake) user@computer:~/openst_adult_demo/spacemake$ ls .
config.yaml project_df.csv  projects  puck_data  species_data
```

(we take the structure from  from the [Adult Mouse tutorial](../examples/adult_mouse/preprocessing_sequencing.md))

From here on, we assume you ran `spacemake`. Please check the [Quick start guide](https://spacemake.readthedocs.io/en/latest/quick-start/index.html#open-st-quick-start) for Open-ST data if you haven't done so.

The typical workflow consists of three steps:

1. **Stitching** of spatial transcriptomics data. This is optional, as it is done automatically by spacemake. Also optionally, tile-scan images can be stitched with the `openst` package.
2. **Merging** of ST and imaging modalities, and pairwise alignment.
3. **Segmentation & DGE creation**, from the merged and pairwise aligned data.

## Spatial (and image) stitching
### Spatial
If you ran `spacemake` with a `puck` with a `coordinate_system`, e.g.:

```yaml
pucks:
  openst:
    coordinate_system: puck_data/openst_coordinate_system.csv
    spot_diameter_um: 0.6
    width_um: 1200
```

**and** a `run_mode` without meshing, e.g.:
```yaml
run_modes:
  openst:
    clean_dge: false
    count_intronic_reads: true
    count_mm_reads: true
    detect_tissue: false
    mesh_data: false
    mesh_type: 'hexagon'
    n_beads: 100000
    polyA_adapter_trimming: true
    spatial_barcode_min_matches: 0.1
```

This means that you will have in your `spacemake` folder a DGE file containing a *spots*-by-*genes* matrix, 
instead of a *bin*-by-*genes* or *cell*-by-*genes* matrices, e.g.:

```sh
(spacemake) user@computer:~/openst_adult_demo/spacemake$ ls projects/openst_demo/processed_data/openst_demo_adult_mouse/illumina/complete_data/dge/dge.*.h5ad
dge.all.polyA_adapter_trimmed.mm_included.spatial_beads_puck_collection.h5ad
```
!!! warning
    The DGE file must **not** contain `mesh_7_hexagon` or similar strings; this means that the spatial spots have been aggregated
    into a regular mesh (in this case, of hexagons with 7 μm side). The files that are necessary for pairwise alignment are
    the barcoded spots (0.6 μm).

Otherwise, if you didn't specify a `coordinate_system` or you only have `mesh`(ed) h5ad DGEs,
you can run this manually, from the `openst` environment (we use `openst_demo` and `openst_demo_adult_mouse`, from the [Adult Mouse tutorial](../examples/adult_mouse/preprocessing_sequencing.md)):

```sh
openst from_spacemake \
     --project-id openst_demo \
     --sample-id openst_demo_adult_mouse \
     spatial_stitch \
     --tile-coordinates fc_2_coordinate_system.csv # replace with the correct coordinate_system
                                                   # this one is available on the Open-ST website at
                                                   # https://rajewsky-lab.github.io/openst/latest/examples/datasets
```

### (Optional) Image
When performing imaging, i.e. as tile-scans, you might need to assemble a single, large image from a collection of smaller tiles (fields-of-view, FOVs). 
We provide code to do this with the microscope setup we used in our manuscript. It is possible to do this using the spacemake file structure. 

For this, please copy the raw data contents into a folder inside the specific sample, e.g.:

```sh
cp images/*.tiff projects/openst_demo/openst_demo/raw_data/images/openst_demo_adult_mouse/.
```

The `destination` folder here (`.../images/openst_demo_adult_mouse`) will contain the image (and metadata) files inside, **these must not be inside nested folders**. 
Otherwise, the `from_spacemake` command will not detect the data. Then, you can run the stitching:

```sh
openst from_spacemake \
     --project-id openst_demo \
     --sample-id openst_demo_adult_mouse \
     image_stitch \
     --imagej-bin ~/Fiji.app/ImageJ-linux64 \ # downloaded externally
     --microscope keyence
```

!!! note
     If you have images that have been stitched outside of the `openst` tools (e.g., automatically by the microscope), you must copy them to:

     ```sh
     cp stitched_image.tiff projects/openst_demo/openst_demo/processed_data/images/Image_Stitched_Composite.tif
     ```

     `openst from_spacemake` assumes the image name is `Image_Stitched_Composite.tif`, therefore the image must be in `tif` format.

## Merging of ST and imaging modalities, pairwise alignment
Once the ST and imaging data have been stitched, you can create a single object containing both modalities:

```sh
openst from_spacemake \
     --project-id openst_demo \
     --sample-id openst_demo_adult_mouse \
     merge_modalities
```

Then, you can run the pairwise alignment, automatically, semiautomatically, or manually. 
You can refer to the [tutorial](pairwise_alignment.md), which provides all details on how to do this. Here we provide some
example commands that replace the ones shown there, to use the `from_spacemake` utility:

=== "Fully automatic alignment"

    ``` sh
    openst from_spacemake \
          --project-id openst_demo \
          --sample-id openst_demo_adult_mouse \
          pairwise_aligner

    # Visualize the alignment
    openst from_spacemake \
          --project-id openst_demo \
          --sample-id openst_demo_adult_mouse \
          manual_pairwise_aligner \
          --spatial-key obsm/spatial_pairwise_aligned_fine \
          --image-key uns/spatial_pairwise_aligned/staining_image_transformed
    ```

=== "Semiautomatic alignment"

    ``` sh
    openst from_spacemake \
          --project-id openst_demo \
          --sample-id openst_demo_adult_mouse \
          pairwise_aligner \
          --only-coarse

    # Visualize the alignment
    # Make sure to save a keypoints.json file
    openst from_spacemake \
          --project-id openst_demo \
          --sample-id openst_demo_adult_mouse \
          manual_pairwise_aligner \
          --spatial-key obsm/spatial_pairwise_aligned_coarse \
          --image-key uns/spatial_pairwise_aligned/staining_image_transformed

    openst from_spacemake \
          --project-id openst_demo \
          --sample-id openst_demo_adult_mouse \
          apply_transform \
          --keypoints-in keypoints.json \
          --spatial-key-in obsm/spatial_pairwise_aligned_coarse \
          --spatial-key-out obsm/spatial_pairwise_aligned_fine \
          --per-tile
    ```

=== "Fully manual alignment"

    ``` sh
    # Visualize the alignment
    # You need to run 'Apply to data' after coarse alignment (key spatial_manual_coarse), then
    # re-render with the transformed data, select keypoints per tile, and
    # make sure to save a keypoints.json file
    openst from_spacemake \
          --project-id openst_demo \
          --sample-id openst_demo_adult_mouse \
          manual_pairwise_aligner \
          --spatial-key obsm/spatial \
          --image-key uns/spatial/staining_image

    openst from_spacemake \
          --project-id openst_demo \
          --sample-id openst_demo_adult_mouse \
          apply_transform \
          --keypoints-in keypoints.json \
          --spatial-key-in obsm/spatial_manual_coarse \
          --spatial-key-out obsm/spatial_manual_fine \
          --per-tile
    ```

Depending on how you ran the alignment, the `spatial-key` with the aligned coordinates might be different.
You can select whatever name you prefer (check the arguments for the commands).


## Segmentation & DGE creation
Once the ST and imaging modalities have been properly aligned, you can segment the images into single cells/nuclei,
and then aggregate the spot locations into individual cells for subsequent analysis. Please adapt the `--image-in` argument
as needed; if you ran fully manual, you can leave as default; if you ran fully automated, you need to specify e.g. `--image-in uns/spatial_pairwise_aligned/staining_image_transformed`.

```sh
openst from_spacemake \
     --project-id openst_demo \
     --sample-id openst_demo_adult_mouse \
     segment \
     --model HE_cellpose_rajewsky # default model for segmebntation of H&E images
```

!!! tip
      You can assess the quality of the segmentation mask using `openst preview`, which leverages `napari` for visualization:

      ```sh
      openst from_spacemake \
            --project-id openst_demo \
            --sample-id openst_demo_adult_mouse \
            preview \
            --image-key uns/spatial/staining_image uns/spatial/staining_image_mask
      ```

      This will create a `napari` window with two image layers. We recommend changing the mask _image_ layer into a
      [_label_ layer](https://napari.org/stable/howtos/layers/labels.html), which is designed for _displaying each integer (ID
      from the segmentation mask) as a different random color, with background rendered as transparent_.

Then, you can create the DGE from the segmentation and ST data. Make sure to populate the `--spatial-key` argument with the
proper key depending on whether you ran (semi)automated or fully manual pairwise aligmment.


```sh
openst from_spacemake \
     --project-id openst_demo \
     --sample-id openst_demo_adult_mouse \
     transcript_assign \
     --spatial-key obsm/spatial_pairwise_aligned_fine # default if you ran (semi)automated pipeline
```

## Expected output
After running all these steps, you will find a `stitched_segmented.h5ad` file inside the `spacemake` folder structure, that you can use for analysis:

```sh
(spacemake) user@computer:~/openst_adult_demo/spacemake$ ls projects/openst_demo/processed_data/openst_demo_adult_mouse/multimodal
stitched_segmented.h5ad
```