# Generating a cell-by-gene matrix
After pairwise alignment, the same coordinate system is shared between the spatial barcodes and the
staining images. 

However, analysis (e.g., clustering, pseudotime, DGE...) is performed on single cells, not on individual capture areas 
(0.6 μm in the current version of the protocol).

Here we show how to aggregate the barcoded spots-by-gene matrix
into a cell-by-genes matrix, where cells are defined from the segmentation mask.

## Segmentation of staining image
To create such a spatial cell-by-gene ($M\times G$) expression matrix, you will first need a segmentation mask.

We efficiently segment cells (or nuclei) from staining images using [cellpose](https://github.com/MouseLand/cellpose).
For the H&E-stained tissue provided in this example, we used our [fine-tuned model](https://github.com/danilexn/openst/blob/main/models/HE_cellpose_rajewsky).
Make sure to download it and save it into a new `models` folder that you need to create under the `openst_e13_demo` main folder.

You can run the segmentation on the previously created `openst_demo_e13_mouse_spatial_beads_puck_collection_aligned.h5ad` file, which
contains the spatial transcriptome coordinates and staining image after coarse+fine pairwise alignment.

```sh
openst segment \
    --adata alignment/openst_demo_e13_mouse_spatial_beads_puck_collection_aligned.h5ad \
    --image-in 'uns/spatial_pairwise_aligned/staining_image_transformed' \
    --output-mask 'uns/spatial_pairwise_aligned/mask_transformed_0px' \
    --model models/HE_cellpose_rajewsky \
    --chunked \
    --gpu \
    --num-workers 8
```

After running this command, the segmentation mask is created and stored in the same `--adata` file, under
the dataset `uns/spatial_pairwise_aligned/mask_transformed_0px`.

## Assigning transcripts to segmented cells
Now, we aggregate the initial barcoded spots-by-gene matrix into a cells-by-genes matrix, by leveraging the
segmentation mask.

!!! note
    In this dataset, we did not specify `--dilate-px`. Since this model segments nuclei from H&E, we are referring
    to nuclei when we talk about cells. We use this terminology for consistency with the rest of the documentation.

This step allows you to aggregate capture spots by segmented cells:

```sh
openst transcript_assign \
    --adata alignment/openst_demo_e13_mouse_spatial_beads_puck_collection_aligned.h5ad \
    --spatial-key spatial_pairwise_aligned_fine \
    --mask-in-adata \
    --mask 'uns/spatial_pairwise_aligned/mask_transformed_0px' \
    --output alignment/openst_demo_e13_mouse_by_cell.h5ad
```

## Expected output
That's it! After running the steps above, you will have that single `alignment/openst_demo_e13_mouse_by_cell.h5ad` file,
which contains the transcriptomic information per segmented cell (nucleus). You can also find [here](https://zenodo),
our output to the one you just generated. 

!!! warning
    Do not expect the files to be exactly the same, as there are several
    steps on the pipeline that are not deterministic (e.g., alignment, segmentation).

Now, the following sections provide some examples of exploratory data analysis that can be performed with this type of data.