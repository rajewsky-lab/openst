# Generating a cell-by-gene matrix
After pairwise alignment, the same coordinate system is shared between the spatial barcodes and the
staining images. 

However, analysis (e.g., clustering, pseudotime, DGE...) is performed on single cells, not on individual capture areas 
(0.6 μm in the current version of the protocol).

So, we show how to aggregate the $N\times G$ matrix ($N$ spots; $G$ genes)
into a $M\times G$ matrix ($M$ segmented cells; $G$ genes), where $N$ maps to $M$ via the segmentation mask.

## Segmentation of staining image
To create such a spatial cell-by-gene ($M\times G$) expression matrix, you will first need a segmentation mask.

We efficiently segment cells (or nuclei) from staining images using [cellpose](https://github.com/MouseLand/cellpose).
We provide a model that we fine-tuned for segmentation of fresh-frozen, H&E-stained tissue,
[here](https://github.com/danilexn/openst/blob/main/models/HE_cellpose_rajewsky).
You can specify any other model that works best for your data -
refer to the [cellpose](https://cellpose.readthedocs.io/en/latest/index.html) documentation.

```sh
openst segment \
    --h5-in spatial_stitched_spots.h5ad \ # after running the pairwise alignment
    --image-in uns/spatial/staining_image \
    --mask-out uns/spatial/staining_image_mask \
    --model HE_cellpose_rajewsky
    # --device cuda \ # uses GPU for segmentation, if available
    # --chunked \ # specify if you run out of GPU memory - segments in chunks
```
By default, segmentation is extended radially 10 pixels. This can be changed with the argument `--dilate-px`.

Make sure to populate the arguments with the values specific to your dataset. Here, we provide `--h5-in` consistent
with the previous steps, `--image-in` and `--mask-out` will read and write the staining and mask inside the Open-ST h5 object,
and `--model` is `HE_cellpose_rajewsky`, the default used in our manuscript. This is the model we recommend for H&E images, and
weights are automatically downloaded. It is also [provided in our repo](https://github.com/rajewsky-lab/openst/blob/main/models/HE_cellpose_rajewsky).
The rest of parameters can be checked with `openst segment --help`.

!!! tip
     **If your sample also contains very large cells** (e.g., adipocytes) that are not segmented with the previous parameters,
     you can perform a second segmentation with a cellpose model, adjusting the diameter parameter.

     ```sh
     openst segment \
          --h5-in spatial_stitched_spots.h5ad \
          --image-in uns/spatial/staining_image \
          --mask-out uns/spatial/staining_image_mask_large \
          --model HE_cellpose_rajewsky \
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

## Assigning transcripts to segmented cells
Now, we aggregate the initial $N\times G$ matrix into an $M\times G$ matrix,
where $N$ maps to $M$ via the segmentation mask.

This step allows you to associate capture spots with segmented cells.

```sh
openst transcript_assign \
    --h5-in spatial_stitched_spots.h5ad \
    --spatial-key obsm/spatial_pairwise_aligned_fine \
    --mask-in uns/spatial/staining_image_mask \
    --h5-out spatial_per_segmented_cell.h5ad
```

In this case, the argument `--mask-in` was set to a single mask, but it can be set to the previously 
introduced _combined_ masks (e.g., `--mask-in uns/spatial/staining_image_mask_combined`).

## Expected output
After running the steps above, you will have a single `h5ad` file, containing the transcriptomic information per segmented cell,
with spatial coordinates compatible with the staining image. The staining image and the segmented image are provided in this object,
so it is possible to visualize it with [squidpy](https://github.com/scverse/squidpy) or [spatialdata](https://github.com/scverse/spatialdata),
among other tools.

So, this concludes the preprocessing of 2D spatial transcriptomics and imaging data
of the Open-ST protocol. Next steps include 3D reconstruction, and
downstream analysis of nD data.