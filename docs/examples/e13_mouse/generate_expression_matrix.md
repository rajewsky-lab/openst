# Segmentation and single-cell quantification
Once the ST and imaging modalities [have been aligned](pairwise_alignment.md), you can segment the images into single cells/nuclei, 
and then aggregate the spot locations into individual cells for subsequent analysis.

## Segmentation of staining image
First, segment the imaging data into single nuclei.
For this dataset, the default `HE_cellpose_rajewsky` works really well.

```sh
openst from_spacemake \
    --project-id openst_demo \
    --sample-id openst_demo_e13_mouse_head \
    segment \
    --image-in 'uns/spatial_pairwise_aligned/staining_image_transformed' \
    --mask-out 'uns/spatial/staining_image_mask' \
    --model models/HE_cellpose_rajewsky \
    --dilate-px 0 \
    --device cuda
```

## Assigning transcripts to segmented cells
Now, aggregate the initial barcoded spots-by-gene matrix into a cells(nuclei)-by-genes matrix, by leveraging the
segmentation mask.

This step allows you to aggregate capture spots by segmented cells (nuclei):

```sh
openst from_spacemake \
    --project-id openst_demo \
    --sample-id openst_demo_e13_mouse_head \
    transcript_assign \
    --spatial-key obsm/spatial_manual_fine
```

## Expected output
That's it! After running the steps above, you will have that single file
which contains the transcriptomic information per segmented cell (nucleus).

!!! warning
    Do not expect the files to be exactly the same to the ones we generated, as there are several
    steps on the pipeline that may not be deterministic (e.g., alignment, segmentation).

Now, the following section provides some examples of exploratory data analysis.