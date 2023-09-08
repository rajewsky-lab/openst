# with experimental spacemake
SMK_DIR='/data/rajewsky/home/dleonpe/projects/bruening_mpi_mm_hypothalamus_latest/data/1_spacemake.bruening_mm_hypo_final_exp/projects/fc_sts_079/processed_data'

ls $SMK_DIR | xargs -n 1 -P 32 -I {} sh -c "python /data/rajewsky/home/dleonpe/projects/openst_paper/code/0_scripts/puck_collection.py \
    --pucks /data/rajewsky/home/dleonpe/projects/bruening_mpi_mm_hypothalamus_latest/data/1_spacemake.bruening_mm_hypo_final_exp/projects/fc_sts_079/processed_data/{}/illumina/complete_data/dge/dge.all.polyA_adapter_trimmed.mm_included.spatial_beads_*.h5ad \
    --puck-coordinates /data/rajewsky/home/dleonpe/projects/bruening_mpi_mm_hypothalamus_latest/notebooks/_misc/align_with_he/novaseq_no_spacing.csv \
    --output /data/rajewsky/home/dleonpe/projects/bruening_mpi_mm_hypothalamus_latest/data/2_downstream.bruening_mm_hypo_final_exp/fc_sts_79/stitched/{}_stitched_spots_for_alignment.h5ad \
    --join-output 'outer' \
    --puck-id-regex '(L[1-4]_tile_[1-2][0-7][0-9][0-9])'" &

SMK_DIR='/data/rajewsky/home/dleonpe/projects/bruening_mpi_mm_hypothalamus_latest/data/1_spacemake.bruening_mm_hypo_final_exp/projects/fc_sts_080/processed_data'
ls $SMK_DIR | xargs -n 1 -P 32 -I {} sh -c "python /data/rajewsky/home/dleonpe/projects/openst_paper/code/0_scripts/puck_collection.py \
    --pucks /data/rajewsky/home/dleonpe/projects/bruening_mpi_mm_hypothalamus_latest/data/1_spacemake.bruening_mm_hypo_final_exp/projects/fc_sts_080/processed_data/{}/illumina/complete_data/dge/dge.all.polyA_adapter_trimmed.mm_included.spatial_beads_*.h5ad \
    --puck-coordinates /data/rajewsky/home/dleonpe/projects/bruening_mpi_mm_hypothalamus_latest/notebooks/_misc/align_with_he/novaseq_no_spacing.csv \
    --output /data/rajewsky/home/dleonpe/projects/bruening_mpi_mm_hypothalamus_latest/data/2_downstream.bruening_mm_hypo_final_exp/fc_sts_80/stitched/{}_stitched_spots_for_alignment.h5ad \
    --join-output 'outer' \
    --puck-id-regex '(L[1-4]_tile_[1-2][0-7][0-9][0-9])'"