#!/bin/bash

# important: run with openst environment (given in environment.yaml)

# # (Optional, microscope dependent) stitch staining image
# openst stitch_image --input-dir $INDIR --output-dir $OUTDIR

# # (Optional) CUT-based preprocessing of the staining image

# # Collect ST tiles into a single file (optional, it is created by spacemake)
# openst create_puck_collection \
#     --pucks dge.all.polyA_adapter_trimmed.mm_included.spatial_beads_*.h5ad \
#     --puck-coordinates /data/rajewsky/home/dleonpe/projects/openst_paper/data/0_misc/novaseq_S4_coordinate_system_FC_009.csv \
#     --output /data/rajewsky/home/dleonpe/projects/openst_paper/data/2_downstream/fc_sts_63/stitched/{}_stitched_spots.h5ad \
#     --join-output 'outer'


# # Segmentation of the (stitched, optionally preprocessed) statining image
# cellpose \
#     --verbose \
#     --use_gpu \
#     --image_path $INDIR/{}/Image_Stitched_Composite.tif \
#     --chan 0 \
#     --pretrained_model $MODEL \
#     --diameter 20 \
#     --flow_threshold 2 \
#     --cellprob_threshold -1 \
#     --save_tif \
#     --savedir $OUTDIR/{} \
#     --no_npy

# Run pairwise alignment
python -m openst pairwise_aligner \
    --image-in /data/rajewsky/home/dleonpe/projects/openst_paper/data/1_images_ln/stitched/fc_sts_63_4/Image_Stitched_Composite.tif \
    --h5-in /data/rajewsky/home/dleonpe/projects/openst_paper/data/2_downstream/fc_sts_63/stitched/fc_sts_063_4_merged_stitched_spots.h5ad \
    --h5-out /data/rajewsky/home/dleonpe/projects/openst_paper/data/2_downstream/fc_sts_63/auto_aligned_he/fc_sts_063_4_stitched_spots_aligned.h5ad \
    --save-image-in-h5 \
    --feature-matcher 'LoFTR' \
    --ransac-coarse-max-trials 1000 \
    --ransac-fine-max-trials 1000 \
    --n-threads 8 \
    --ransac-coarse-residual-threshold 2 \
    --ransac-coarse-min-samples 2 \
    --metadata-out /data/rajewsky/home/dleonpe/projects/openst_paper/data/2_downstream/fc_sts_63/auto_aligned_he/fc_sts_063_4_alignment_metadata.pkl \
    --threshold-counts-coarse 2

# # Assign transcripts to cells using the segmentation mask. Needs metadata alignment
# openst assign_transcripts_to_mask