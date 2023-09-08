# with experimental spacemake

python /data/rajewsky/home/dleonpe/projects/openst_paper/repos/openst/src/alignment/pairwise_aligner.py \
    --image-in /data/rajewsky/home/dleonpe/projects/openst_paper/data/1_images_ln/stitched/fc_sts_63_4/Image_Stitched_Composite.tif \
    --h5-in /data/rajewsky/home/dleonpe/projects/openst_paper/data/2_downstream/fc_sts_63/stitched/fc_sts_063_4_merged_stitched_spots.h5ad \
    --h5-out /data/rajewsky/home/dleonpe/projects/openst_paper/data/2_downstream/fc_sts_63/auto_aligned_he/fc_sts_063_4_stitched_spots_aligned_loftr.h5ad \
    --save-image-in-h5 \
    --feature-matcher 'LoFTR' \
    --ransac-coarse-max-trials 1000 \
    --ransac-fine-max-trials 1000 \
    --n-threads 8 \
    --ransac-coarse-residual-threshold 2 \
    --ransac-coarse-min-samples 2 \
    --metadata-out /data/rajewsky/home/dleonpe/projects/openst_paper/data/2_downstream/fc_sts_63/auto_aligned_he/fc_sts_063_4_alignment_metadata_loftr.pkl \
    --threshold-counts-coarse 2