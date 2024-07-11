# API of Open-ST tools

Here we show the description of all commands available in the `openst` tools, which can be run as:

```bash
openst <subcommand>
```

## `barcode_preprocessing`
Convert spatial barcode raw data into tabular files with barcodes and spatial coordinates.

Usage:
```text
openst barcode_preprocessing [-h] --fastq-in FASTQ_IN --tilecoords-out TILECOORDS_OUT --out-suffix OUT_SUFFIX [--out-prefix OUT_PREFIX] [--crop-seq CROP_SEQ] [--rev-comp] [--single-tile] [--unsorted]

options:
  -h, --help            show this help message and exit
  --fastq-in FASTQ_IN   Path to the fastq file
  --tilecoords-out TILECOORDS_OUT
                        Directory where output files will be written to
  --out-suffix OUT_SUFFIX
                        Suffix added to the name of the output files (i.e., extension)
  --out-prefix OUT_PREFIX
                        (Optional) Prefix added to the name of the output files. Default: ""
  --crop-seq CROP_SEQ   (Optional) A 'python-style' slice, used to crop input sequences. Default: ":"
  --rev-comp            (Optional) Apply reverse complementary after sequence cropping
  --single-tile         (Optional) set if it is guarranteed that the input .fastq(.gz) file contains only a tile
  --unsorted            (Optional) set when file is unsorted respect to tiles; might be slower
```

## `image_stitch`
Stitch image fields of view into a single image. Currently, it only supports `--microscope keyence`, for the
default microscopy setup used in our paper.

Usage:
```text
openst image_stitch [-h] --image-indir IMAGE_INDIR --image-out IMAGE_OUT --imagej-bin IMAGEJ_BIN --microscope {keyence} [--no-run] [--rerun] [--metadata METADATA] [--join-zstack-regex JOIN_ZSTACK_REGEX]

options:
  -h, --help            show this help message and exit
  --image-indir IMAGE_INDIR
                        path to collection of images
  --image-out IMAGE_OUT
                        path where to save the image (must be a filename)
  --imagej-bin IMAGEJ_BIN
                        path to the ImageJ/Fiji executable. Must have the Grid Collection plugin available!
  --microscope {keyence}
                        microscope model or imaging strategy that was used for imaging
  --no-run              If set, do not run ImageJ, but return the command line
  --rerun               If set, runs stitching even when the output file exists
  --metadata METADATA   Path where the metadata will be stored. If not specified, metadata is not saved. Default: ""
  --join-zstack-regex JOIN_ZSTACK_REGEX
                        When non empty, this specifies how to find the Z location from the individual filename and will create a z-stack from single images. Example regex: 'Image_([0-9]*)_Z([0-9]*)_CH1.tif'. Default: ""
```

## `image_preprocess`
Restoration of imaging data with CUT model, as in Open-ST paper.

Usage:
```text
openst image_preprocess [-h] [--h5-in H5_IN] [--image-in IMAGE_IN] [--image-out IMAGE_OUT] [--tile-size-px TILE_SIZE_PX] [--model MODEL] [--device {cpu,cuda}] [--num-workers NUM_WORKERS]

options:
  -h, --help            show this help message and exit
  --h5-in H5_IN         If set, image is loaded from the Open-ST h5 object (key in --image-in), and retored image is saved there (to the key --image-out). Default: ""
  --image-in IMAGE_IN   Key or path to the input image. Default: ""
  --image-out IMAGE_OUT
                        Key or path where the restored image will be written into. Default: ""
  --tile-size-px TILE_SIZE_PX
                        The input image is split into squared tiles of side `--tile-size-px`, for inference.Larger values avoid boundary effects, but require more memory. Default: 512
  --model MODEL         CUT model used for image restoration. Default: "HE_CUT_rajewsky"
  --device {cpu,cuda}   Device used to run CUT restoration model. Can be ['cpu', 'cuda']. Default: "cpu"
  --num-workers NUM_WORKERS
                        Number of CPU workers for parallel processing. Default: -1
```

## `spatial_stitch`
Stitch Open-ST h5 tile objects into a single Open-ST h5 object.

Usage:
```text
openst spatial_stitch [-h] --tiles TILES [TILES ...] --tile-coordinates TILE_COORDINATES [--tile-id TILE_ID [TILE_ID ...]] [--tile-id-regex TILE_ID_REGEX] --h5-out H5_OUT [--tile-id-key TILE_ID_KEY]
                             [--merge-output {same,unique,first,only}] [--join-output {inner,outer}] [--no-reset-index] [--no-transform] [--metadata METADATA]

options:
  -h, --help            show this help message and exit
  --tiles TILES [TILES ...]
                        Path to h5 file, one per tile - separated by space
  --tile-coordinates TILE_COORDINATES
                        Path to the coordinate system file
  --tile-id TILE_ID [TILE_ID ...]
                        (Mandatory if --tile-id-regex is not specified) Per tile file specified in --tiles, each entry in --tile-id maps a tile file to the tile IDs under the first column of the --tile-
                        coordinates file.
  --tile-id-regex TILE_ID_REGEX
                        (Mandatory if --tile-id is not specified) "Regex to find tile id in file names, instead of specifying a list in --tile-id. Default: "(L[1-4][a-b]_tile_[1-2][0-7][0-9][0-9])"
  --h5-out H5_OUT       Where the stitched spatial object will be written to
  --tile-id-key TILE_ID_KEY
                        Key of the h5 file (under /obs) where tile IDs are stored. If != 'tile_id', a new categorical column of this name will be generated for consistency. Default: "tile_id"
  --merge-output {same,unique,first,only}
                        how to merge tiles, can be "same", "unique", "first", "only". Default: "same"
  --join-output {inner,outer}
                        how to join tiles, can be "inner", "outer". Default: "outer"
  --no-reset-index      If set, do not reset the obs_name index of the combined spatial object as 'obs_name:<tile_id_key>'; keep original 'obs_name'.
  --no-transform        If set, spatial coordinates are not transformed - just combine tiles into a single spatial object
  --metadata METADATA   (Optional) Path where the metadata will be stored. If not specified, metadata is not saved. Default: ""
```

## `merge_modalities`
Merges spatial locations (as points) and images of Open-ST data into a single h5 object.

Usage
```text
openst merge_modalities [-h] --h5-in H5_IN --image-in IMAGE_IN [--image-key IMAGE_KEY]

options:
  -h, --help            show this help message and exit
  --h5-in H5_IN         Input Open-ST h5 object
  --image-in IMAGE_IN   Image that will be loaded and written into the Open-ST h5 object
  --image-key IMAGE_KEY
                        Key in the Open-ST h5 object where the image will be saved. Default: "uns/spatial/staining_image"
```

## `pairwise_aligner`
Automatic pairwise alignment of transcript locations to imaging data.

Usage:
```text
openst pairwise_aligner [-h] --h5-in H5_IN [--image-in IMAGE_IN] [--metadata METADATA] [--only-coarse] [--rescale-factor-coarse RESCALE_FACTOR_COARSE] [--threshold-counts-coarse THRESHOLD_COUNTS_COARSE]
                               [--pseudoimage-size-coarse PSEUDOIMAGE_SIZE_COARSE] [--ransac-coarse-min-samples RANSAC_COARSE_MIN_SAMPLES] [--ransac-coarse-residual-threshold RANSAC_COARSE_RESIDUAL_THRESHOLD]
                               [--ransac-coarse-max-trials RANSAC_COARSE_MAX_TRIALS] [--genes-coarse GENES_COARSE [GENES_COARSE ...]] [--rescale-factor-fine RESCALE_FACTOR_FINE]
                               [--tissue-masking-gaussian-sigma TISSUE_MASKING_GAUSSIAN_SIGMA] [--fine-registration-gaussian-sigma FINE_REGISTRATION_GAUSSIAN_SIGMA]
                               [--threshold-counts-fine THRESHOLD_COUNTS_FINE] [--pseudoimage-size-fine PSEUDOIMAGE_SIZE_FINE] [--ransac-fine-min-samples RANSAC_FINE_MIN_SAMPLES]
                               [--ransac-fine-residual-threshold RANSAC_FINE_RESIDUAL_THRESHOLD] [--ransac-fine-max-trials RANSAC_FINE_MAX_TRIALS] [--fine-min-matches FINE_MIN_MATCHES]
                               [--genes-fine GENES_FINE [GENES_FINE ...]] [--mask-tissue] [--keep-black-background] [--feature-matcher {LoFTR,SIFT,KeyNet}] [--fiducial-model FIDUCIAL_MODEL]
                               [--num-workers NUM_WORKERS] [--device {cpu,cuda}]

options:
  -h, --help            show this help message and exit

Data (required):
  --h5-in H5_IN         Path to the merged Open-ST h5 object containing spatial coordinates and images

Data (optional):
  --image-in IMAGE_IN   Key to the image used as the 'destination' during pairwise alignment. Default: "uns/spatial/staining_image"
  --metadata METADATA   Path where the metadata will be stored. If not specified, metadata is not saved. Default: ""

Coarse registration parameters:
  --only-coarse         If selected, only the coarse alignment stage will run
  --rescale-factor-coarse RESCALE_FACTOR_COARSE
                        Rescaling factor for the input image (1:factor), used during coarse pairwise alignment. Default: 20
  --threshold-counts-coarse THRESHOLD_COUNTS_COARSE
                        Only spatial coordinates with counts larger than this number will be kept for pseudoimage rendering during coarse alignment. Default: 1
  --pseudoimage-size-coarse PSEUDOIMAGE_SIZE_COARSE
                        Size (in pixels) of the pseudoimage during coarse alignment. Default: 500
  --ransac-coarse-min-samples RANSAC_COARSE_MIN_SAMPLES
                        'min_samples' parameter of RANSAC, during coarse registration. Default: 3
  --ransac-coarse-residual-threshold RANSAC_COARSE_RESIDUAL_THRESHOLD
                        'residual_threshold' parameter of RANSAC, during coarse registration. Default: 2
  --ransac-coarse-max-trials RANSAC_COARSE_MAX_TRIALS
                        Times RANSAC will run (x1000 iterations) during coarse registration. Default: 2
  --genes-coarse GENES_COARSE [GENES_COARSE ...]
                        Genes used for plotting the pseudoimage during the coarse alignment phase. Default: None

Fine registration parameters:
  --rescale-factor-fine RESCALE_FACTOR_FINE
                        Rescaling factor for the input image (1:factor), used during fine pairwise alignment. Default: 10
  --tissue-masking-gaussian-sigma TISSUE_MASKING_GAUSSIAN_SIGMA
                        The gaussian blur sigma used during the isolation of the tissue on the HE (preprocessing). Default: 5
  --fine-registration-gaussian-sigma FINE_REGISTRATION_GAUSSIAN_SIGMA
                        Gaussian blur used on all modalities during fine registration. Default: 2
  --threshold-counts-fine THRESHOLD_COUNTS_FINE
                        Only spatial coordinates with counts larger than this number will be kept for pseudoimage rendering during fine alignment. Default: 0
  --pseudoimage-size-fine PSEUDOIMAGE_SIZE_FINE
                        Size (in pixels) of the pseudoimage during fine alignment. Default: 2000
  --ransac-fine-min-samples RANSAC_FINE_MIN_SAMPLES
                        'min_samples' parameter of RANSAC, during fine registration. Default: 3
  --ransac-fine-residual-threshold RANSAC_FINE_RESIDUAL_THRESHOLD
                        'residual_threshold' parameter of RANSAC, during fine registration. Default: 2
  --ransac-fine-max-trials RANSAC_FINE_MAX_TRIALS
                        Times RANSAC will run (x1000 iterations) during fine registration. Default: 2
  --fine-min-matches FINE_MIN_MATCHES
                        Minimum number of matching keypoints between modalities during fine alignment. Default: 50
  --genes-fine GENES_FINE [GENES_FINE ...]
                        Genes used for plotting the pseudoimage during the fine alignment phase. Default: None

Image preprocessing parameters:
  --mask-tissue         Tissue (imaging modality) is masked from the background for the feature detection
  --keep-black-background
                        Whether to set the background of the imaging modalities to white, after tissue masking

Feature model parameters:
  --feature-matcher {LoFTR,SIFT,KeyNet}
                        Feature matching algorithm. Default: "LoFTR"
  --fiducial-model FIDUCIAL_MODEL
                        Path to a object detection model (YOLO) to detect fiducial markers. Default: ""

Computational parameters:
  --num-workers NUM_WORKERS
                        Number of CPU workers for parallel processing. Default: 1
  --device {cpu,cuda}   Device used to run feature matching model. Can be ['cpu', 'cuda']. Default: "cpu"
```

## `apply_transform`
Apply a precomputed transformation matrix to the specified coordinates of an Open-ST h5 object

Usage:
```text
openst apply_transform [-h] --keypoints-in KEYPOINTS_IN --h5-in H5_IN [--per-tile] [--spatial-key-in SPATIAL_KEY_IN] [--spatial-key-out SPATIAL_KEY_OUT]

options:
  -h, --help            show this help message and exit
  --keypoints-in KEYPOINTS_IN
                        Path to the json file containing keypoints.
  --h5-in H5_IN         Path to the input h5ad file containing spatial coordinates
  --per-tile            (Optional) If set, transformations are applied per tile, from their keypoints. Otherwise, a single transform is computed for all tiles.
  --spatial-key-in SPATIAL_KEY_IN
                        Key of the Open-ST h5 object where the input spatial coordinates are read from. Default: "obsm/spatial_pairwise_aligned_coarse"
  --spatial-key-out SPATIAL_KEY_OUT
                        Key of the Open-ST h5 object where the transformed spatial coordinates are written into. Default: "obsm/spatial_pairwise_aligned_fine"
```

## `manual_pairwise_aligner`
GUI for manual alignment of Open-ST data

```text
openst manual_pairwise_aligner [-h] [--h5-in H5_IN] [--spatial-key SPATIAL_KEY] [--image-key IMAGE_KEY]

options:
  -h, --help            show this help message and exit
  --h5-in H5_IN         Path to the input h5ad file containing spatial coordinates. Default: ""
  --spatial-key SPATIAL_KEY
                        Path in the h5ad file to the spatial coordinates. Default: ""
  --image-key IMAGE_KEY
                        Path in the h5ad file to the image. Default: ""
```
## `segment`
Image (or pseudoimage)-based segmentation with cellpose and (optional) radial extension


Usage:
```text
openst segment [-h] [--image-in IMAGE_IN] [--h5-in H5_IN] --mask-out MASK_OUT [--rna-segment] [--model MODEL] [--flow-threshold FLOW_THRESHOLD] [--cellprob-threshold CELLPROB_THRESHOLD]
                      [--diameter DIAMETER] [--chunk-size CHUNK_SIZE] [--chunked] [--max-image-pixels MAX_IMAGE_PIXELS] [--device {cpu,cuda}] [--dilate-px DILATE_PX] [--outline-px OUTLINE_PX] [--mask-tissue]
                      [--tissue-masking-gaussian-sigma TISSUE_MASKING_GAUSSIAN_SIGMA] [--keep-black-background] [--rna-segment-spatial-coord-key RNA_SEGMENT_SPATIAL_COORD_KEY]
                      [--rna-segment-input-resolution RNA_SEGMENT_INPUT_RESOLUTION] [--rna-segment-render-scale RNA_SEGMENT_RENDER_SCALE] [--rna-segment-render-sigma RNA_SEGMENT_RENDER_SIGMA]
                      [--rna-segment-output-resolution RNA_SEGMENT_OUTPUT_RESOLUTION] [--num-workers NUM_WORKERS] [--metadata METADATA]

options:
  -h, --help            show this help message and exit
  --image-in IMAGE_IN   Key in the Open-ST h5 object (when --h5-in is specified) or path to the file where the mask will be loaded from
  --h5-in H5_IN         If specified, image is loaded from h5 (from key --image-in). Segmentation mask is saved there (to --mask-out). Default: ""
  --mask-out MASK_OUT   Key in the Open-ST h5 object (when --h5-in is specified) or path to the file where the mask will be written into
  --rna-segment         Performs segmentation based on local RNA density pseudoimages from sequencing data, instead of using a staining image. This assumes coordinates in microns (can be transformed with
                        --rna-segment-input-resolution)
  --model MODEL         cellpose model - either a path or a valid string to pretrained model. Default: ""
  --flow-threshold FLOW_THRESHOLD
                        cellpose's 'flow_threshold' parameter. Default: 0.5
  --cellprob-threshold CELLPROB_THRESHOLD
                        cellpose's 'cellprob_threshold' parameter. Default: 0
  --diameter DIAMETER   cellpose's 'diameter' parameter. Default: 20
  --chunk-size CHUNK_SIZE
                        When prediction of the mask runs in separate chunks, this is the chunk square size (in pixels). Default: 512
  --chunked             If set, segmentation is computed at non-overlapping chunks of size '--chunk-size'
  --max-image-pixels MAX_IMAGE_PIXELS
                        Upper bound for number of pixels in the images (prevents exception when opening very large images). Default: 933120000
  --device {cpu,cuda}   Device used to run the segmentation model. Can be ['cpu', 'cuda']. Default: "cpu"
  --dilate-px DILATE_PX
                        Pixels the outlines of the segmentation mask will be extended. Default: 10
  --outline-px OUTLINE_PX
                        Objects will be represented as px-width outlines (only if >0). Default: 0
  --mask-tissue         Tissue (imaging modality) is masked from the background before segmentation.
  --tissue-masking-gaussian-sigma TISSUE_MASKING_GAUSSIAN_SIGMA
                        The gaussian blur sigma used during the isolation of the tissue on the staining image. Default: 5
  --keep-black-background
                        Whether to set the background of the imaging modalities to white after tissue masking
  --rna-segment-spatial-coord-key RNA_SEGMENT_SPATIAL_COORD_KEY
                        Path to the spatial coordinates inside the spatial object (e.g., 'obsm/spatial'). Default: "obsm/spatial"
  --rna-segment-input-resolution RNA_SEGMENT_INPUT_RESOLUTION
                        Spatial resolution of the input coordinates (retrieved from --rna-segment-spatial-coord-key). If it is in microns, leave as 1. If it is in pixels, specify the pixel to micron conversion
                        factor. Default: 1
  --rna-segment-render-scale RNA_SEGMENT_RENDER_SCALE
                        Size of bins for computing the binning (in microns). For Open-ST v1, we recommend a value of 2. Default: 2
  --rna-segment-render-sigma RNA_SEGMENT_RENDER_SIGMA
                        Smoothing factor applied to the RNA pseudoimage (higher values lead to smoother images). Default: 1
  --rna-segment-output-resolution RNA_SEGMENT_OUTPUT_RESOLUTION
                        Final resolution (micron/pixel) for the segmentation mask. Default: 0.6
  --num-workers NUM_WORKERS
                        Number of CPU workers when --chunked is specified. Default: -1
  --metadata METADATA   Path where the metadata will be stored. If not specified, metadata is not saved. Warning: a report (via openst report) cannot be generated without metadata! Default: ""
```

## `segment_merge`
Merge two segmentation masks into one

Usage:
```text
openst segment_merge [-h] --h5-in H5_IN --mask-in MASK_IN MASK_IN --mask-out MASK_OUT [--chunk-size CHUNK_SIZE] [--chunked] [--num-workers NUM_WORKERS]

options:
  -h, --help            show this help message and exit
  --h5-in H5_IN         If set, masks are loaded from the Open-ST h5 object (key in --mask-in), and segmentation is saved there (to the key under --mask-out). Default: ""
  --mask-in MASK_IN MASK_IN
                        Path to the input segmentation masks - two of them!
  --mask-out MASK_OUT   Path (file or h5) where the merged mask will be saved
  --chunk-size CHUNK_SIZE
                        When prediction of the mask runs in separate chunks, this is the chunk square size (in pixels). Default: 512
  --chunked             If set, segmentation is computed at non-overlapping chunks of size '--chunk-size'
  --num-workers NUM_WORKERS
                        Number of CPU workers when --chunked is specified. Default: -1
```

## `transcript_assign`
Aggregate transcripts into segmented cells.

Usage:

```text
openst transcript_assign [-h] --h5-in H5_IN --mask-in MASK_IN --spatial-key SPATIAL_KEY --h5-out H5_OUT [--mask-from-file] [--max-image-pixels MAX_IMAGE_PIXELS] [--shuffle-umi] [--metadata METADATA]

options:
  -h, --help            show this help message and exit
  --h5-in H5_IN         Path to an already aligned Open-ST h5 object
  --mask-in MASK_IN     Path to image mask - a key in the Open-ST h5 object. Or, can be an image stored separately in the filesystem (when --mask-from-file is specified) Image data and ST coordinates must be
                        pairwise aligned (implicit for the case of RNA-based segmentation)
  --spatial-key SPATIAL_KEY
                        Key in the Open-ST h5 object where the aligned coordinates are stored, e.g. 'spatial_pairwise_aligned_coarse' (after using 'openst pairwise_aligner')
  --h5-out H5_OUT       Path where the segmented Open-ST h5 object will be written into
  --mask-from-file      If set, the image mask is loaded from an external file
  --max-image-pixels MAX_IMAGE_PIXELS
                        Upper bound for number of pixels in the images (prevents exception when opening very large images). Default: 933120000
  --shuffle-umi         If set, UMI locations will be shuffled. This can be used as a baseline for feature selection.
  --metadata METADATA   Path where the metadata will be stored. If not specified, metadata is not saved. Warning: a report (via openst report) cannot be generated without metadata! Default: ""
```

## `pseudoimage`
Generate pseudoimages of Open-ST RNA data and visualize using napari.

Usage:
```text
openst pseudoimage [-h] --h5-in H5_IN [--spatial-coord-key SPATIAL_COORD_KEY] [--input-resolution INPUT_RESOLUTION] [--render-scale RENDER_SCALE] [--render-sigma RENDER_SIGMA]
                          [--output-resolution OUTPUT_RESOLUTION]

options:
  -h, --help            show this help message and exit
  --h5-in H5_IN         Necessary to create the pseudoimage
  --spatial-coord-key SPATIAL_COORD_KEY
                        Path to the spatial coordinates inside the spatial object (e.g., 'obsm/spatial')
  --input-resolution INPUT_RESOLUTION
                        Spatial resolution of the input coordinates (retrieved from --spatial-coord-key). If it is in microns, leave as 1. If it is in pixels, specify the pixel to micron conversion factor. Default: 1
  --render-scale RENDER_SCALE
                        Size of bins for computing the binning (in microns). For Open-ST v1, we recommend a value of 2. Default: 2
  --render-sigma RENDER_SIGMA
                        Smoothing factor applied to the RNA pseudoimage (higher values lead to smoother images). Default: 1
  --output-resolution OUTPUT_RESOLUTION
                        Final resolution (micron/pixel) for the pseudoimage. Default: 0.6
```

## `preview`
Preview locations (as points) and images of Open-ST data.

Usage:
```text
openst preview [-h] --h5-in H5_IN [--file-structure] [--spatial-coord-keys SPATIAL_COORD_KEYS [SPATIAL_COORD_KEYS ...]] [--image-keys IMAGE_KEYS [IMAGE_KEYS ...]]
                      [--pseudoimage-keys PSEUDOIMAGE_KEYS [PSEUDOIMAGE_KEYS ...]] [--spatial-coord-resampling SPATIAL_COORD_RESAMPLING [SPATIAL_COORD_RESAMPLING ...]]
                      [--image-resampling IMAGE_RESAMPLING [IMAGE_RESAMPLING ...]] [--pseudoimage-units-to-um PSEUDOIMAGE_UNITS_TO_UM [PSEUDOIMAGE_UNITS_TO_UM ...]]

options:
  -h, --help            show this help message and exit
  --h5-in H5_IN         Necessary to create the pseudoimage
  --file-structure      If set, will not open a visualization screen but will return the tree structure of the h5 file
  --spatial-coord-keys SPATIAL_COORD_KEYS [SPATIAL_COORD_KEYS ...]
                        Path to the spatial coordinates inside the spatial object (e.g., 'obsm/spatial'). Can be one or many (separated by space)
  --image-keys IMAGE_KEYS [IMAGE_KEYS ...]
                        Path to the image to be visualized. Can be one or many (separated by space)
  --pseudoimage-keys PSEUDOIMAGE_KEYS [PSEUDOIMAGE_KEYS ...]
                        Path to the spatial coordinates inside the spatial object to visualize as pseudoimage. Can be one or many (separated by space)
  --spatial-coord-resampling SPATIAL_COORD_RESAMPLING [SPATIAL_COORD_RESAMPLING ...]
                        Will load every n-th point. Can be one (same for all spatial-coords) or many (1-to-1 mapping to the spatial-coord list). Default: [1]
  --image-resampling IMAGE_RESAMPLING [IMAGE_RESAMPLING ...]
                        Will load every n-th pixel. Can be one (same for all images) or many (1-to-1 mapping to the image list). Default: [1]
  --pseudoimage-units-to-um PSEUDOIMAGE_UNITS_TO_UM [PSEUDOIMAGE_UNITS_TO_UM ...]
                        Conversion factor from spatial units to micron, before rendering the pseudoimage. Can be one (same for all images) or many (1-to-1 mapping to the image list). Default: [1.0]
```

## `report`
Generate HTML reports from metadata files (json).

Usage:
```text
openst report [-h] --metadata METADATA --html-out HTML_OUT

options:
  -h, --help           show this help message and exit
  --metadata METADATA  Path to the metadata file (json)
  --html-out HTML_OUT  Path where the output HTML file will be created
```

## `from_spacemake`
Run openst commands using spacemake file structure. You need to specify one `subcommand` from above, with the respective arguments.

Usage:

```text
openst from_spacemake [-h] --project-id PROJECT_ID --sample-id SAMPLE_ID [--run-mode RUN_MODE] subcommand [params]

options:
  -h, --help            show this help message and exit
  --project-id PROJECT_ID
                        From spacemake's project_df, this is the project_id string
  --sample-id SAMPLE_ID
                        From spacemake's project_df, this is the sample_id string
  --run-mode RUN_MODE   When a sample has multiple run_mode(s), you must specify one
```

This command populates the arguments for the subcommands automatically. In the table below, you can find
which of this arguments are populated automatically. `...` indicates that the value will change
depending on the `--project-id` and `--sample-id` configuration.

| subcommand | populated arguments |
| ---- | ---- |
| `spatial_stitch` | `--h5-out ... --tiles ... --tile-id ...` |
| `image_stitch` | `--image-indir ... --image-out ...` |
| `segment_merge` | `--h5-in ... --mask-out uns/spatial/staining_image_mask_merged` |
| `segment` | `--h5-in ... --image-in uns/spatial/staining_image --mask-out uns/spatial/staining_image_mask`|
| `transcript_assign` | `--h5-in ... --mask-in uns/spatial/staining_image_mask --h5-out ...` |
| `merge_modalities` | `--h5-in ... --image-in ...` |
| `manual_pairwise_aligner` | `--h5-in ...` |
| `apply_transform` | `--h5-in ... `|
| `pseudoimage` | `--h5-in ... `|
| `preview` | `--h5-in ... `|
| `pairwise_aligner` | `--h5-in ... `|

!!! note
    You can override any of these values by providing them explicitly.