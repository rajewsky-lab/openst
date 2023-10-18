# 3D reconstruction from serial sections of spatial transcriptomics and H&E images
In this section, we will guide you through the process of creating a 3D reconstruction from serial sections of spatial transcriptomics (ST) and H&E images using the Spatial Transcriptomics ImgLib2/Imaging Project (STIM, v0.2.0). This reconstruction allows you to gain a comprehensive understanding of your biological samples in three dimensions.

## Creation of csv files
Use the provided script

```bash
openst to_3d_registration --args --metadata=<where_to_write_metadata>
```

## Conversion to n5 format
Convert the coordinate and gene expression information of these datasets into the n5 format, which is optimized for efficient image processing, using the st-resave function.
```bash
STIMBINS="/home/dleonpe/data/bin"
STIMINFILES="/data/rajewsky/home/dleonpe/projects/openst_paper/data/2_downstream/fc_sts_63/aligned_sections/1_input"
STIMOUTFILES="/data/rajewsky/home/dleonpe/projects/openst_paper/data/2_downstream/fc_sts_63/aligned_sections/2_stim_dataset"
FNAME="stitched_spots_merged_aligned_10px_GAN_segmented.h5ad."

$STIMBINS/st-resave \
    -i "$STIMINFILES/fc_sts_63_2_${FNAME}locations.csv,$STIMINFILES/fc_sts_63_2_${FNAME}genes.csv,fc_sts_63_02" \
    -i "$STIMINFILES/fc_sts_63_3_${FNAME}locations.csv,$STIMINFILES/fc_sts_63_3_${FNAME}genes.csv,fc_sts_63_03" \
    -o "$STIMOUTFILES/fc_sts_63.n5" \
    --normalize
```

## Pairwise alignment
Utilize the st-align-pairs function to perform pairwise alignment of three sections below and above each section (r=3). This function creates image channels of gene expression for prespecified genes, aggregated per cell as a Gauss rendering around centroids, parametrized with a smoothness factor.

```bash
STIMBINS="/home/dleonpe/data/bin"
STIMOUTFILES="/data/rajewsky/home/dleonpe/projects/openst_paper/data/2_downstream/fc_sts_63/aligned_sections/2_stim_dataset"

# 2. Run pairwise alignment
$STIMBINS/st-align-pairs \
     -i "$STIMOUTFILES/fc_sts_63.n5" \
     --scale 0.03 \
     -r 3 \
     --hidePairwiseRendering \
     --overwrite \
     -sf 4.0 \
     --minNumInliers 15 \
     --numGenes 0 \
     -g 'KRT6A,KRT6B,S100A2,LYZ,CD74,IGKC,IGHG1,IGHA1,JCHAIN,CD74,AMTN'
     #--numGenes 20 \

```

## Feature Filtering and Global Alignment
Filter the resulting set of feature matches between pairs of sections using an affine model. Configure the st-align-pairs function with appropriate parameters, such as --minNumInliers 15, --scale 0.03, and -sf 4.0 (smoothness factor).
```bash
STIMBINS="/home/dleonpe/data/bin"
STIMOUTFILES="/data/rajewsky/home/dleonpe/projects/openst_paper/data/2_downstream/fc_sts_63/aligned_sections/2_stim_dataset"

$STIMBINS/st-align-global \
     -i "$STIMOUTFILES/fc_sts_63.n5" \
     --skipICP \
     -g 'KRT6A,KRT6B,S100A2,LYZ,CD74,IGHG1,IGHA1,JCHAIN,CD74,AMTN'
```

## Conversion to h5ad Format
Convert the n5 container back to the h5ad format for subsequent downstream analyses. This will also transfer the transformation models from the ST alignment onto the preprocessed and background-removed H&E images. This script will output the aligned spatial coordinates and an image volume that can be used for subsequent 3D visualization, of spatial transcriptomics and H&E staining in a common coordinate system.
```bash
openst from_3d_registration --args --metadata=<where_to_write_metadata>
```

You can quickly generate a HTML report to visualize the alignment quality. This will provide, for instance, a volumetric rendering on your browser using the channels (genes) selected for registration. This will also visualize the sections individually, to assess the deformations per section after registration, as well as the deformation fields.

```bash
openst report --metadata=<where_to_write_metadata> --output=<path_to_html_file>
```

With these steps completed, you will have successfully reconstructed a 3D representation of your biological samples, integrating spatial transcriptomics data and H&E images. This 3D reconstruction provides valuable insights into the spatial distribution of gene expression within your samples and enhances your understanding of complex biological structures.

## 3D visualization with ParaView
(explain)