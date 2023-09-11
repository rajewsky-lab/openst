# Open-ST: Computational Protocol

## Introduction

(explain the aim)

(explain the general steps)

1. Preprocessing of spatial transcriptomics data
2. Preprocessing of imaging data
3. Alignment of imaging and spatial transcriptomics
4. Assigning transcripts to segmented cells
5. (Optional) 3D reconstruction from serial sections of spatial transcriptomics and H&E images

## Requirements

What to download, what to install

## 1. Preprocessing sequencing data

(we use spacemake from the fastq reads)

### 1.1 Running spacemake
You can find more instructions on the spacemake documentation. Here we provide an example, refer there for the complete documentation.

### 1.2 Quality control

We check the sample quality (e.g., PCR bias), and the stitching files

## 2. Preprocessing of imaging data
We'll start with the initial processing of H&E-stained tissue sections, as a prerequisite to aligning spatial transcriptomics data with tissue staining images.

### 2.1 Stitching (microscope-dependent)
To begin, stitch together the tile-scan images of your H&E-stained tissue sections using the Grid/Collection stitching plugin included in Fiji 1.53t. This will create a composite image of the entire section.

*Warning: this step depends on the microscope used for imaging. In our protocol, we use the (XXX microscope), and provide code to perform the stitching on any computer. Refer to the documentation or vendor of your microscope for the stitching of tile-scans*

```bash
# placeholder for code
```

### 2.2 (Optional) Style transfer
Optionally, you can perform style transfer on your tile-scan images using a custom Contrastive Unpaired Translation (CUT) model. This can help equalize the style between sections and remove artifacts.

```bash
# placeholder for code
```

### 2.3 Segmentation of staining image
Next, segment the nuclei in your images using Cellpose 2.2. We provide a model optimized for segmentation of fresh-frozen, H&E-stained tissue. You can specify any other model that works best for your data; refer to the cellpose documentation.

```bash
# placeholder for code
```

#### Extra: segmentation of very large cells
If your samples contain very large cells that cannot be segmented with the provided H&E model (e.g., adipocytes), you can perform a second round of segmentation with a cellpose model, adjusting the diameter parameter.

```bash
# placeholder for code
```
Then, you can extend and combine the segmentation masks of both diameter configurations.

```bash
# placeholder for code
```

## 3. Alignment of imaging and spatial transcriptomics
In order to assign transcripts to the nuclei segmented from the staining images, the pairwise alignment between the imaging and spatial transcriptomics modality must be performed. For this, we provide software that allows to:
1. Creation of pseudoimages from the spatial transcriptomics data.
2. Two-step (coarse, then fine) alignment of H&E images to pseudoimages of ST data.
3. Manual curation and refinement of the alignment for more precision (~0.5 µm error), by leveraging fiducial markers visible from both modalities.

### 3.1 Pairwise alignment
What to run
#### 3.2 Quality Control
With the pweaver report
#### 3.3 (Optional) Refinement of alignment
Use the provided notebook that spawns a napari environment for visual assessment and manual refinement of the alignment, if necessary.
```bash
# placeholder for code
```

## 4. Assigning transcripts to segmented cells
Finally, create spatial cell-by-gene expression matrices by aggregating the initial NxG matrix into an MxG matrix, where N maps to M via the segmentation mask. This step allows you to associate capture spots with segmented cells.

```bash
openst run --verbose 1
```

That concludes the preprocessing of imaging data in the Open-ST computational protocol. 

## 5. (Optional) 3D reconstruction from serial sections of spatial transcriptomics and H&E images
In this section, we will guide you through the process of creating a 3D reconstruction from serial sections of spatial transcriptomics (ST) and H&E images using the Spatial Transcriptomics ImgLib2/Imaging Project (STIM, v0.2.0). This reconstruction allows you to gain a comprehensive understanding of your biological samples in three dimensions.

### 5.1 Creation of csv files
Use the provided script

### 5.2 Conversion to n5 format
Convert the coordinate and gene expression information of these datasets into the n5 format, which is optimized for efficient image processing, using the st-resave function.

### 5.3 Pairwise alignment
Utilize the st-align-pairs function to perform pairwise alignment of three sections below and above each section (r=3). This function creates image channels of gene expression for prespecified genes, aggregated per cell as a Gauss rendering around centroids, parametrized with a smoothness factor.

```bash
# Placeholder for pairwise alignment command
```

### 5.4 Feature Filtering and Global Alignment
Filter the resulting set of feature matches between pairs of sections using an affine model. Configure the st-align-pairs function with appropriate parameters, such as --minNumInliers 15, --scale 0.03, and -sf 4.0 (smoothness factor).
```bash
# Placeholder for pairwise alignment command
```

### 5.5 Conversion to h5ad Format
Convert the n5 container back to the h5ad format for subsequent downstream analyses.
```bash
# Placeholder for pairwise alignment command
```

### 5.6 Alignment of H&E images
Transfer the resultant affine transformation matrices obtained from the ST alignment onto the preprocessed and background-removed H&E images. These images can be rescaled to any 1:X factor. This script will also create an aligned imaging volume that can be used for subsequent 3D visualization.

```bash
# Placeholder for pairwise alignment command
```

With these steps completed, you will have successfully reconstructed a 3D representation of your biological samples, integrating spatial transcriptomics data and H&E images. This 3D reconstruction provides valuable insights into the spatial distribution of gene expression within your samples and enhances your understanding of complex biological structures.

### 5.7 3D visualization with ParaView
(explain)

## FAQ