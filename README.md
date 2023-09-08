# Automatic Pairwise Alignment of Spatial Transcriptomics and Imaging Data (open-ST)
Author: Daniel León-Periñán @ N.Rajewsky Lab (BIMSB)
Date: August 30, 2023

## Description
This script provides a solution for performing automatic pairwise alignment between spatial transcriptomics (STS) data and imaging data; particularly, for our open-ST method. The primary purpose of this alignment is to overlay spatial transcriptomics data onto imaging data, enabling spatial comparison and analysis.

## Usage Example

```bash

python pairwise_aligner.py --image-in input_image.jpg --h5-in input_data.h5ad --h5-out aligned_data.h5ad
                     --metadata-out metadata.pkl --save-image-in-h5 --rescale-factor-coarse 20
                     --rescale-factor-fine 5 --tissue-masking-gaussian-sigma 5 --fine-registration-gaussian-sigma 2
                     --set-white-background --threshold-counts-coarse 1 --threshold-counts-fine 0
                     --pseudoimage-size-coarse 4000 --pseudoimage-size-fine 6000 --ransac-coarse-min-samples 3
                     --ransac-coarse-residual-threshold 2 --ransac-coarse-max-trials 10000 --ransac-fine-min-samples 10
                     --ransac-fine-residual-threshold 2 --ransac-fine-max-trials 10000 --max-image-pixels 933120000
                     --feature-matcher 'SIFT' --n-threads 2
```

### Features
- Performs automatic pairwise alignment between spatial transcriptomics (STS) data and imaging data.
- Aligns STS data onto imaging data for spatial comparison and analysis.
- Supports various alignment parameters such as rescaling factors, Gaussian blur settings, and more.
- Allows the saving of aligned images and metadata into the output h5ad file.
- Provides options for reporting and logging the alignment process.

### Requirements
- Python 3.10
- Third-party libraries: scikit-image, numpy, h5py, torch, scikit-learn, kornia, anndata

### Installation
- Clone or download this repository.
  ```bash
  git clone https://github.com/rajewsky-lab/openst
  ```
- Install the provided environment using conda
  ```bash
  conda create -n environment.yaml
  ```

### License
This project is licensed under the MIT License - see the LICENSE file for details.