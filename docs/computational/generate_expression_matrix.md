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
    --h5-in <path_to_aligned_h5ad> \
    --image-in <image_in_path> \
    --mask-out <mask_out_path> \
    --model <path>/HE_cellpose_rajewsky
    # --device cuda \ # uses GPU for segmentation, if available
    # --chunked \ # specify if you run out of GPU memory - segments in chunks
```
By default, segmentation is extended radially 10 pixels. This can be changed with the argument `--dilate-px`.

Make sure to replace the placeholders (`<...>`). For instance,
`<path_to_aligned_h5ad>` is the full path to the `h5ad` file [after pairwise alignment](pairwise_alignment.md#expected-output); 
`<image_in_path>` is the path to the image - a path to a file, or a location inside the `h5ad` file,
like `'uns/spatial_pairwise_aligned/staining_image_transformed'` (*our recommendation*).
`<mask_out_path>` is the location where the segmentation mask will be saved - can be a file or a location in the `h5ad` file,
like `uns/spatial_pairwise_aligned/mask_transformed_10px` (*our recommendation*). The `<model_path>` for the parameter `--model`
is the name or location of the cellpose model weights.

We recommend using the [model provided in our repo](https://github.com/rajewsky-lab/openst/blob/main/models/HE_cellpose_rajewsky)
for segmentation of H&E images The rest of parameters can be checked with `openst segment --help`.

!!! tip
     **If your sample also contains very large cells** (e.g., adipocytes) that are not segmented with the previous parameters,
     you can perform a second segmentation with a cellpose model, adjusting the diameter parameter.

     ```sh
     openst segment \
          --h5-in <path_to_aligned_h5ad> \
          --image-in <image_in_path> \
          --mask-out <mask_out_path_larger> \
          --model <path>/HE_cellpose_rajewsky \
          --dilate-px 50 \
          --diameter 50 # diameter for the larger cell type
     ```

     Replace the placeholders (`<...>`) as before; in this case, the placeholder `<mask_out_path_larger>` must be different from the
     `<mask_out_path>` provided above.
     
     And then, you can combine the segmentation masks of both diameter configurations.
     This command will apply an "AND" between all images, to only preserve mask of non-overlapping,
     with the hierarchy provided in the `--image-in` argument (first has higher priority).

     ```sh
     openst segment_merge \
          --h5-in <path_to_aligned_h5ad> \
          --mask-in <mask_a> <mask_b>
          --mask-out <mask_combined>
     ```

     Replace the placeholders (`<...>`) as before; in this case, the placeholder `<mask_a>`, `<mask_b>`... must correspond
     to the placeholders `<mask_out_path>`, `<mask_out_path_larger>`...

## Assigning transcripts to segmented cells
Now, we aggregate the initial $N\times G$ matrix into an $M\times G$ matrix,
where $N$ maps to $M$ via the segmentation mask.

This step allows you to associate capture spots with segmented cells.

```sh
openst transcript_assign \
    --h5-in <path_to_aligned_h5ad> \
    --spatial-key spatial_pairwise_aligned_fine \
    --mask-in<mask_out_path> \
    --h5-out <path_to_sc_h5ad>
```

Replace the placeholders (`<...>`) as before; in this case, the placeholder `<mask_in_path>` must be set to
be equal to the `<mask_out_path>` (or `<mask_combined>` if you ran multiple segmentation); also, `<path_to_sc_h5ad>`
must be set to a valid path and filename where the output cell-by-gene matrix (not barcode-by-cell) will be written.

## Expected output
After running the steps above, you will have a single `h5ad` file, containing the transcriptomic information per segmented cell,
with spatial coordinates compatible with the staining image. The staining image and the segmented image are provided in this object,
so it is possible to visualize it with [squidpy](https://github.com/scverse/squidpy) or [spatialdata](https://github.com/scverse/spatialdata),
among other tools.

So, this concludes the preprocessing of 2D spatial transcriptomics and imaging data
of the Open-ST protocol. Next steps include 3D reconstruction, and
downstream analysis of nD data.