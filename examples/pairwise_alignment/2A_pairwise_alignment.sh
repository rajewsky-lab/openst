# with experimental spacemake
SMK_DIR='/data/rajewsky/home/dleonpe/projects/bruening_mpi_mm_hypothalamus_latest/data/1_spacemake.bruening_mm_hypo_final_exp/projects/fc_sts_079/processed_data'

ls $SMK_DIR | xargs -n 1 -P 16 -I {} sh -c "python /data/rajewsky/home/dleonpe/projects/openst_paper/repos/openst/src/alignment/pairwise_aligner.py \
    --image-in /data/rajewsky/home/dleonpe/projects/bruening_mpi_mm_hypothalamus_latest/data/1_images.bruening_mm_hypo_final_exp/stitched/{}/Image_Stitched_Composite.tif \
    --h5-in /data/rajewsky/home/dleonpe/projects/bruening_mpi_mm_hypothalamus_latest/data/2_downstream.bruening_mm_hypo_final_exp/fc_sts_79/stitched/{}_stitched_spots_for_alignment.h5ad \
    --h5-out /data/rajewsky/home/dleonpe/projects/bruening_mpi_mm_hypothalamus_latest/data/2_downstream.bruening_mm_hypo_final_exp/fc_sts_79/auto_aligned_he/{}_stitched_spots_aligned_loftr.h5ad \
    --save-image-in-h5 \
    --feature-matcher 'LoFTR' \
    --ransac-coarse-max-trials 1000 \
    --ransac-fine-max-trials 1000 \
    --n-threads 8 \
    --ransac-coarse-residual-threshold 2 \
    --ransac-coarse-min-samples 2 \
    --metadata-out /data/rajewsky/home/dleonpe/projects/bruening_mpi_mm_hypothalamus_latest/data/2_downstream.bruening_mm_hypo_final_exp/fc_sts_80/auto_aligned_he/{}_alignment_metadata_loftr.pkl \
    --threshold-counts-coarse 2" &

SMK_DIR='/data/rajewsky/home/dleonpe/projects/bruening_mpi_mm_hypothalamus_latest/data/1_spacemake.bruening_mm_hypo_final_exp/projects/fc_sts_080/processed_data'
ls $SMK_DIR | xargs -n 1 -P 16 -I {} sh -c "python /data/rajewsky/home/dleonpe/projects/openst_paper/repos/openst/src/alignment/pairwise_aligner.py \
    --image-in /data/rajewsky/home/dleonpe/projects/bruening_mpi_mm_hypothalamus_latest/data/1_images.bruening_mm_hypo_final_exp/stitched/{}/Image_Stitched_Composite.tif \
    --h5-in /data/rajewsky/home/dleonpe/projects/bruening_mpi_mm_hypothalamus_latest/data/2_downstream.bruening_mm_hypo_final_exp/fc_sts_80/stitched/{}_stitched_spots_for_alignment.h5ad \
    --h5-out /data/rajewsky/home/dleonpe/projects/bruening_mpi_mm_hypothalamus_latest/data/2_downstream.bruening_mm_hypo_final_exp/fc_sts_80/auto_aligned_he/{}_stitched_spots_aligned_loftr.h5ad \
    --save-image-in-h5 \
    --feature-matcher 'LoFTR' \
    --ransac-coarse-max-trials 1000 \
    --ransac-fine-max-trials 1000 \
    --n-threads 8 \
    --ransac-coarse-residual-threshold 2 \
    --ransac-coarse-min-samples 2 \
    --metadata-out /data/rajewsky/home/dleonpe/projects/bruening_mpi_mm_hypothalamus_latest/data/2_downstream.bruening_mm_hypo_final_exp/fc_sts_80/auto_aligned_he/{}_alignment_metadata_loftr.pkl \
    --threshold-counts-coarse 2"

# python /data/rajewsky/home/dleonpe/projects/openst_paper/repos/openst/src/alignment/pairwise_aligner.py \
#     --image-in /data/rajewsky/home/dleonpe/projects/bruening_mpi_mm_hypothalamus_latest/data/1_images.bruening_mm_hypo_final_exp/stitched/fc_sts_079_8/Image_Stitched_Composite.tif \
#     --h5-in /data/rajewsky/home/dleonpe/projects/bruening_mpi_mm_hypothalamus_latest/data/2_downstream.bruening_mm_hypo_final_exp/fc_sts_79/stitched/fc_sts_079_8_stitched_spots_for_alignment.h5ad \
#     --h5-out /data/rajewsky/home/dleonpe/projects/bruening_mpi_mm_hypothalamus_latest/data/2_downstream.bruening_mm_hypo_final_exp/fc_sts_79/auto_aligned_he/fc_sts_079_8_stitched_spots_aligned.h5ad \
#     --save-image-in-h5 \
#     --feature-matcher 'SIFT' \
#     --ransac-coarse-max-trials 1000 \
#     --ransac-fine-max-trials 1000 \
#     --n-threads 32 \
#     --only-coarse \
#     --ransac-coarse-residual-threshold 2 \
#     --ransac-coarse-min-samples 2 \
#     --metadata-out /data/rajewsky/home/dleonpe/projects/bruening_mpi_mm_hypothalamus_latest/data/2_downstream.bruening_mm_hypo_final_exp/fc_sts_79/auto_aligned_he/fc_sts_079_8_alignment_metadata.pkl \
#     --threshold-counts-coarse 2