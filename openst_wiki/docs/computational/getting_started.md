# Getting started

(explain the aim)

(explain the general steps)

1. Preprocessing of spatial transcriptomics data
2. Preprocessing of imaging data
3. Alignment of imaging and spatial transcriptomics
4. Assigning transcripts to segmented cells
5. (Optional) 3D reconstruction from serial sections of spatial transcriptomics and H&E images

## Installation
The processing of open-ST data relies on several computational steps that can be installed using *conda* or *mamba* (our recommendation). Make sure that this is installed on your machine, following the adequate [instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) for your operating system. 

Once you have installed *conda* or *mamba*, create an environment using the provided [environment.yaml](https://github.com/danilexn/openst/blob/main/environment.yaml) as a template. To do this, download this file, browse to the folder where it was downloaded, and create the environment with the following command.

```bash
conda ... (from yaml)
```

At any time, you can activate the openst environment via:

```bash
conda activate openst
```

Then, another critical dependency for the preprocessing of the sequencing data is spacemake. We recommend following the [official instructions](https://spacemake.readthedocs.io/en/latest/). Since spacemake and openst processing are respectively sequential, we recommend creating separate environments for spacemake and openst.

## Optional dependencies
There are two dependencies that need to be installed. On the one hand, [Fiji](https://imagej.net/software/fiji/downloads) might be necessary for the stitching of tile-scan, tissue staining images (see ). On the other hand, in the case of performing 3D registration of consecutive serial sections, we recommend [STIM](https://github.com/PreibischLab/STIM).