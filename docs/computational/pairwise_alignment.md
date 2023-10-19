# Align image to transcriptome
In the previous step, the transcriptomic reads were processed and mapped in tissue space with spacemake.
Now, in order to assign transcripts to cells from the staining images, we perform a pairwise alignment 
between the imaging and spatial transcriptomics modality.

`openst` provides tools to automatically or manually carry out this pairwise alignment operation. For this,
our tool generates *pseudoimages* of the spatial transcriptomics data, such that computer vision algorithms
can be used for the alignment of true staining images and the ST modalities.

The alignment workflow consists of two steps, that can be performed [automatically] or [manually]:

1. Coarse alignment of H&E images to pseudoimages of ST data - at low resolution.
2. Fine alignment, using fiducial marks detected at both modalities for very precise alignment.

[automatically]: #automated-workflow
[manually]: #manual-workflow

## Required input data
For automatic and manual alignment, two inputs are required: (1) a stitched tile-scan of 
the staining image (see [Preprocessing of imaging](computational/preprocessing_imaging.md)), and (2) a
single [h5ad] file containing all the [barcoded tiles] of a sample 
(see [Preprocessing of sequencing](computational/preprocessing_sequencing.md)). 

!!! warning
    Remember that spacemake generates one file per [barcoded tile], 
    and it will be necessary to perform [stitching of tiles] to
    obtain a single file with all tiles for a sample.

### Image modality
The expected input for the image modality is a single tiff file of the whole tile-scan.
Other formats, i.e., `jpeg` and `png` are also compatible. Make sure that the image is the full-resolution
and not a downsampled version. As a rule of thumb, at least a few fiducial markers must be clearly distinguishable.

=== "Good resolution"

    put image

=== "Bad resolution"

    put image

### Spatial transcriptomics modality
The expected input for the spatial transcriptomics modality is a single [h5ad] file containing all
the [barcoded tiles] of a sample. 

!!! warning
    Remember that spacemake generates one file per [barcoded tile], and it will be necessary to perform
    [stitching of tiles] to obtain a single file with all tiles for a sample.

The structure of this single file must follow the [h5ad] standard. In particular, we expect that there is a
column under `obs` called `tile_id`, telling which tile a barcode belongs to. Also, we expect a column 
`total_counts` indicating the total count of unique transcripts per barcode. Finally, we expect that there is
a barcodes-by-2d matrix under `obsm` (with the name `spatial`), containing the spatial coordinates of the barcoded
spots.

[h5ad]: https://anndata.readthedocs.io/en/latest/fileformat-prose.html

## Automated workflow
If you want to save time (😉), we provide a script that performs coarse and fine steps of alignment 
automatically, by leveraging computer vision algorithms. To do so, make sure that you have the [necessary
input data](#required-input-data); then, open a termina, type and run the following command (just an example):

```bash
openst pairwise_aligner \
    --image-in <path_image>/Image_Stitched_Composite.tif \
    --h5-in <path>/<id>_stitched_spots.h5ad \
    --h5-out <path>/<id>_stitched_spots_aligned.h5ad \
    --metadata-out <path>/<id>_alignment_metadata.json \
    --save-image-in-h5
```

Make sure to replace the placeholders (`<...>`). For instance, `<path_image>` in the `--image-in` command 
with the folder containing the [stitched images](computational/preprocessing_imaging.md) of the current sample, 
`<path>` in the `--h5-in` and `--h5-out` arguments  to contain the folder containing the 
[stitched barcoded tiles](computational/preprocessing_sequencing.md), and `<id>` with the `sample_id` as 
defined in the spacemake project. **Importantly**, make sure to specify a path where the metadata output file 
should be created via `--metadata-out`; this will be useful for a visual assessment of whether
automated alignment worked or not.

If you want to run only the coarse phase of the pairwise alignment (i.e., to run the fine
alignment [yourself](#manual-workflow)), you can specify the argument `--only-coarse`. If you have a CUDA-compatible
GPU in your machine, you can specify the argument `--device cuda` to accelerate feature detection and matching. 

!!! note
    For aligning STS to H&E-stained tissues, **we recommend** leaving the arguments with the default values. 
    Anyway, you can get a full list of configurable parameters by running `openst pairwise_aligner --help`.

### Visual assessment of alignment
Before proceeding to generating a cell-by-gene matrix, **we strongly recommend** visually assessing the alignment.
For this, make sure that you specified a file name in the `--metadata-out`; then, you can open a terminal and run
the following command:

```sh
openst report --metadata=<metadata_file> --output=<path_to_html_file>
```

Make sure to replace the placeholders (`<...>`) with the path to the metadata file (under `--metadata` argument), and
the desired path and filename of the HTML report that will be generated (under `--output` argument). Running this command
will create a HTML report containing images of the STS and staining image before and after alignment (coarse and/or fine,
depending on the configuration).

## Manual workflow
If the automatic alignment is not successful, or you prefer to run the alignment yourself, we provide a
Graphical User Interface (GUI) tool that allows to perform manual alignment. This GUI only needs an input file
generated by the `openst` pipeline containing the spatial transcriptome (barcoded spots-by-genes) and the staining
image, for visual reference. Such a file is generated upon running `openst pairwise_aligner` (if you'd like to
refine an automatic alignment). Alternatively, you can start from the [same files as in the automated workflow](#required-input-data),
and merge them with the command:

```sh
openst manual_pairwise_aligner 
openst manual_pairwise_aligner \
    --prepare-data \
    --h5-in <path>/<id>_stitched_spots_aligned.h5ad \
    --image-in <path_image>/Image_Stitched_Composite.tif
```

Make sure to replace the placeholders (`<...>`). For instance, `<path_image>` in the `--image-in` command 
with the folder containing the [stitched images](computational/preprocessing_imaging.md) of the current sample, 
`<path>` in the `--h5-in` argument to contain the folder containing the 
[stitched barcoded tiles](computational/preprocessing_sequencing.md), and `<id>` with the `sample_id` as 
defined in the spacemake project. 

!!! warning
    When `--h5-out` is not provided, the image will be stored as a new layer in the
    file specified at `--h5-in`.

Once you have prepared the input file, you can run the GUI with the following command (no arguments are required)

```sh
openst manual_pairwise_aligner_gui
```

This GUI can be used to align the images from scratch, or to validate and refine the results from
automatic coarse and/or fine alignment. In any case, the GUI is the same: we provide a video walthrough
on how to use this GUI tool for all these cases. In summary, the GUI will allow you to create a `json` file
containing a list of pairs of corresponding points between the staining and spatial transcriptome coordinates.
This can be used later to transform the coordinates of the spatial trancriptomics to match the coordinate system of the image.

---

:fontawesome-brands-youtube:{ style="color: #EE0F0F" }
__[Walkthrough of the GUI for manual alignment]__ by @danilexn – :octicons-clock-24:
5m – Learn how to align STS and imaging data in a step-by-step guide.

  [Walkthrough of the GUI for manual alignment]: https://www.youtube.com

---

Once the `json` list of point correspondences has been generated with the GUI, you can run the following command to transform
the coordinates of the spatial transcriptomics to match the coordinates of the staining image. We provide three alternative commands,
depending on whether you chose to align from scratch, make a manual fine alignment from a (manual or automated) coarse alignment, 
or refine an automated (or manual) fine alignment

=== "Coarse from raw"

    ``` sh
    openst manual_pairwise_aligner \
        --keypoints-json <path_to_keypoints.json> \
        --h5-in <path_to_sts.h5ad> \
        --coarse
    ```

=== "Fine from coarse"

    ``` sh
    openst manual_pairwise_aligner \
        --keypoints-json <path_to_keypoints.json> \
        --h5-in <path_to_sts.h5ad> \
        --fine
    ```

=== "Refine from fine"

    ``` sh
    openst manual_pairwise_aligner \
        --keypoints-json <path_to_keypoints.json> \
        --h5-in <path_to_sts.h5ad> \
        --refine
    ```

Make sure to replace the placeholders (`<...>`). For instance,
`<path_to_keypoints.json>` is the `json` file generated with the GUI, and  in the `--h5-in` argument to contain the folder containing the 
[stitched barcoded tiles](computational/preprocessing_sequencing.md), and `<path_to_sts.h5ad>` is the path to the h5ad file that was loaded
with the GUI. Running the command above will generate a new [obsm] layer containing the transformed spatial transcriptome coordinates

[obsm]: https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.obsm.html

## Expected output
After running the automatic or manual alignment, you must have a single `h5ad` file, containing the transformed spatial coordinates.
This will be used in the following step to aggregate the transcripts by a spatially-corresponding cell, in order to get a cell-by-gene
matrix that can be used in later downstream analysis.