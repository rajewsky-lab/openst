# Segmentation and single-cell quantification
Once the ST and imaging modalities [have been aligned](pairwise_alignment.md#option-a-automated-alignment), 
you can segment the images into single cells/nuclei, and then aggregate the spot locations into individual cells 
for subsequent analysis.

## Cell segmentation from tissue image
Let's create a new Open-ST h5 object containing a cell-by-gene expression matrix. First, you will need a cell 
(or nuclear) segmentation mask.

```sh
openst from_spacemake \
     --project-id openst_demo_project \
     --sample-id openst_demo_sample \
     segment \
     --model HE_cellpose_rajewsky # default model for segmentation of H&E images
```

=== "From (semi)automatic alignment"

     ```sh
     openst from_spacemake \
          --project-id openst_demo_project \
          --sample-id openst_demo_sample \
          segment \
          --model HE_cellpose_rajewsky \
          --image-in uns/spatial_pairwise_aligned/staining_image_transformed
     ```

=== "From manual alignment"

     ```sh
     openst from_spacemake \
          --project-id openst_demo_project \
          --sample-id openst_demo_sample \
          segment \
          --model HE_cellpose_rajewsky
     ```

We segment cells (or nuclei) from staining images using [cellpose](https://github.com/MouseLand/cellpose).
We provide a model that we fine-tuned for segmentation of fresh-frozen, H&E-stained tissue,
[here](http://bimsbstatic.mdc-berlin.de/rajewsky/openst-public-data/models/HE_cellpose_rajewsky), but you can use
any other model (e.g., pretrained from cellpose, like `cyto2` or `nuclei`, or your own).
Also, by default, segmentation is extended radially 10 pixels (see `--dilate-px`), to account for cytoplasm surrounding
the nucleus as a first approximation of cell shape (might hold or not depending on the tissue).

??? question "I want to segment very small and very large cells..."

     You can perform an additional round of segmentation by, e.g., adjusting the diameter parameter.

     ```sh
     openst from_spacemake \
          --project-id openst_demo_project \
          --sample-id openst_demo_sample \
          segment \
          --mask-out uns/spatial/staining_image_mask_large \
          --dilate-px 50 \
          --diameter 50 # diameter for the larger cell type
     ```

     In this case, we changed `--mask-out` to a different key, so we can keep both masks inside the Open-ST h5 object.
     
     Then, you can combine the segmentation masks of both diameter configurations.
     This command will apply an "AND" between all images, to only preserve mask of non-overlapping,
     with the hierarchy provided in the `--image-in` argument (first has higher priority).

     ```sh
     openst segment_merge \
          --h5-in spatial_stitched_spots.h5ad \
          --mask-in uns/spatial/staining_image_mask uns/spatial/staining_image_mask_large
          --mask-out uns/spatial/staining_image_mask_combined
     ```

## Quality control of segmentation
You can assess the quality of segmentation with `openst preview`:

=== "From (semi)automatic alignment"

     ```sh
     openst from_spacemake \
          --project-id openst_demo_project \
          --sample-id openst_demo_sample \
          preview \
          --image-key uns/spatial_pairwise_aligned/staining_image_transformed uns/spatial/staining_image_mask
     ```

=== "From manual alignment"

     ```sh
     openst from_spacemake \
          --project-id openst_demo_project \
          --sample-id openst_demo_sample \
          preview \
          --image-key uns/spatial/staining_image uns/spatial/staining_image_mask
     ```

This will create a `napari` window with two image layers. Change the mask _image_ layer into a
[_label_ layer](https://napari.org/stable/howtos/layers/labels.html), which is designed for _displaying each integer (ID
from the segmentation mask) as a different random color, with background rendered as transparent_.

If you are satistied with the quality of the segmentation, **you are all set to continue with single-cell quantification**.

## Single-cell quantification

Then, you can create a single file containing the transcriptomic information aggregated into (segmented) single-cells.

=== "From automatic alignment"

     ```sh
     openst from_spacemake \
          --project-id openst_demo_project \
          --sample-id openst_demo_sample \
          transcript_assign \
          --spatial-key obsm/spatial_pairwise_aligned_fine \
          --mask-in uns/spatial_pairwise_aligned/staining_image_transformed
     ```

=== "From semiautomatic alignment"

     ```sh
     openst from_spacemake \
          --project-id openst_demo_project \
          --sample-id openst_demo_sample \
          transcript_assign \
          --spatial-key obsm/spatial_manual_fine \
          --mask-in uns/spatial_pairwise_aligned/staining_image_transformed
     ```

=== "From manual alignment"

     ```sh
     openst from_spacemake \
          --project-id openst_demo_project \
          --sample-id openst_demo_sample \
          transcript_assign \
          --spatial-key obsm/spatial_manual_fine
     ```

## Expected output
After the steps above, you will have a single `h5ad` file with transcriptomic information per segmented cell,
with spatial coordinates aligned to the staining image. The staining image and the segmented image are provided in this object,
so it is possible to visualize it with [squidpy](https://github.com/scverse/squidpy) or [spatialdata](https://github.com/scverse/spatialdata),
among other tools.

!!! warning
     In the Open-ST h5 object, the cell with ID 0 will correspond to the background. Please remove it before
     proceeding with analysis.

This concludes the preprocessing of 2D spatial transcriptomics and imaging data
of the Open-ST protocol. Next steps include 3D reconstruction, and
downstream analysis of nD data.