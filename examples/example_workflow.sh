#!/bin/bash

# author: Daniel Leon-Perinan (Rajewsky Lab, 2023)
# this file uses openst to preprocess STS and imaging data, 
# yields a single-cell h5ad file (UMI aggregated by segmentation mask)

# placeholders:
# <path> = path where results from this script will be saved
# <path_image> = path to the folder containing the corresponding staining image
# <path_spacemake_dge> = spacemake-generated path to the dge files of this sample
# <id> = the sample id, as in spacemake

### STEP 0 ###
# activate the openst environment
# e.g., conda activate openst

### STEP 1 ###
# stitching of the individual tiles with a coordinate system
openst spatial_stitch \
    --tiles <path_spacemake_dge>/dge.all.polyA_adapter_trimmed.mm_included.spatial_beads_*.h5ad \
    --tile-coordinates <path>/novaseq_S4_coordinate_system_FC_010.csv \
    --output <path>/<id>_stitched_spots.h5ad \
    --join-output 'outer' \
    --tile-id-regex '(L[1-4]_tile_[1-2][0-7][0-9][0-9])'

### STEP 2 ###
# pairwise alignment of the coordinates
# we only run coarse alignment
# this will integrate the staining image into the file (nice for segmentation)
openst pairwise_aligner \
    --image-in <path_image>/Image_Stitched_Composite.tif \
    --h5-in <path>/<id>_stitched_spots.h5ad \
    --h5-out <path>/<id>_stitched_spots_aligned.h5ad \
    --save-image-in-h5 \
    --feature-matcher 'LoFTR' \
    --ransac-coarse-max-trials 1000 \
    --ransac-fine-max-trials 1000 \
    --n-threads 8 \
    --ransac-coarse-residual-threshold 2 \
    --ransac-coarse-min-samples 2 \
    --only-coarse \
    --device cuda \
    --metadata-out <path>/<id>_alignment_metadata.json \
    --threshold-counts-coarse 2

### STEP 3 ###
# requires manual input from the user!
# we run the refinement with the gui
# if a display is not set from ssh session, run ssh -X <host>, then e.g.,
# export DISPLAY=localhost:10.0
openst manual_pairwise_aligner_gui

### STEP 4 ###
# transfer the refined keypoints to the fine aligned data
openst manual_pairwise_aligner \
    --keypoints-json <path>/<id>_keypoints.json \
    --h5-in <path>/<id>_stitched_spots_aligned.h5ad \
    --fine

### STEP 5 ###
# calculate a segmentation mask for the aligned staining image
# optionally, we apply a dilation of 10px (increases UMIs, but might increase background signal)
openst segment \
    --adata <path>/<id>_stitched_spots_aligned.h5ad \
    --image-in uns/spatial_pairwise_aligned/staining_image_transformed \
    --output-mask uns/spatial_pairwise_aligned/mask_transformed_10px \
    --model <path>/HE_cellpose_rajewsky \
    --cellprob-threshold -1 \
    --flow-threshold 2 \
    --chunked \
    --gpu \
    --num-workers 32 \
    --dilate-px 10

### STEP 6 ###
# assign the individual barcoded spots (0.5 micron) into the segmentation mask
openst transcript_assign \
    --adata <path>/<id>_stitched_spots_aligned.h5ad \
    --spatial-key spatial_pairwise_aligned_fine \
    --mask-in-adata \
    --mask uns/spatial_pairwise_aligned/mask_transformed_10px \
    --output <path>/<id>_segmented_cells_fine.h5ad