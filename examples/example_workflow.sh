# 0. cd into a root spacemake folder (i.e., containing project_df.csv and config.yaml)
# cd ...

# 1. stitching of the spatial dge
openst from_spacemake --project-id fc_sts_76 --sample-id fc_sts_76_3 spatial_stitch --tile-coordinates fc_2_coordinate_system.csv

# 2. copy images into raw directory

# 3. stitching of the images
# TODO: write tutorial for how to download ImageJ
openst from_spacemake --project-id fc_sts_76 --sample-id fc_sts_76_3 image_stitch --imagej-bin ~/Fiji.app/ImageJ-linux64 --microscope keyence

# 4. merge modalities
openst from_spacemake --project-id fc_sts_76 --sample-id fc_sts_76_3 merge_modalities

# 5.A pairwise alignment, automatic
openst from_spacemake --project-id fc_sts_76 --sample-id fc_sts_76_3 pairwise_aligner --device cuda --mask-tissue --only-coarse

# 5.B pairwise alignment, manual
# we load the GUI, then select points for coarse and fine registration. Do not forget to save the keypoints.json file!
openst from_spacemake --project-id fc_sts_76 --sample-id fc_sts_76_3 manual_pairwise_aligner --spatial-key obsm/spatial --image-key uns/spatial/staining_image
# we apply the transform to the data, per tile (cannot be done in GUI)
openst from_spacemake --project-id fc_sts_76 --sample-id fc_sts_76_3 apply_transform --keypoints keypoints.json --spatial-key-in obsm/spatial_pairwise_aligned_coarse --spatial-key-out obsm/spatial_pairwise_aligned_fine --per-tile

# 6. segmentation
openst from_spacemake --project-id fc_sts_76 --sample-id fc_sts_76_3 segment --device cuda --model HE_cellpose_rajewsky

# 7. create dge
openst from_spacemake --project-id fc_sts_76 --sample-id fc_sts_76_3 transcript_assign --spatial-key obsm/spatial_pairwise_aligned_fine
