# Preprocessing of imaging data
If you used our same experimental setup for imaging, you will obtain individual images
as a result of performing imaging of the stained tissue in a tile-scan fashion. From this point,
the Open-ST pipeline expects to have a single image for the whole tile-scan, rather than individual files
per tile.

This step depends on the microscope used for imaging. In [our implementation](../experimental/library_preparation.md#he-staining-and-imaging), we used a
Keyence BZ-X710 inverted fluorescence phase contrast microscope, for which we provide open-source
image stitching code - runs in any computer, independent of the software provided by the manufacturer.

!!! warning
     If you use a different microscope, please refer to the documentation of your microscope
     for how to stitch tile-scans into a single image. The Open-ST pipeline expects one single file in either
     `tiff`, `jpeg` or `png` formats.

## Stitching with the script
When using our same setup, we provide a script that automatically handles the imaging data generated 
by the Keyence microscope and leverages the Grid/Collection stitching plugin included in
[Fiji](https://imagej.net/software/fiji/downloads) to create a single, composite image of the tile-scan.
Open a terminal, and run the following command:

```bash
openst image_stitch \
     --microscope='keyence' \
     --imagej-bin=<path_to_fiji_or_imagej> \
     --image-indir=<path_to_tiles> \
     --image-out=Image_Stitched_Composite.tif # just an example
```
For instance, `<path_to_fiji_or_imagej>` is the path where the [Fiji](https://imagej.net/software/fiji/downloads) executable is;
`<path_to_tiles>` is the folder where the microscope has saved the individual image files with the tile scan
and the `.bcl` file with metadata (in the case of _keyence_ microscope).

!!! question
     If you don't know how to specify the `<path_to_fiji_or_imagej>`, please follow the official instructions provided
     for [Running Headless](https://imagej.net/learn/headless). For example, under linux you can download it using

     ```bash
     cd ~
     wget "https://downloads.imagej.net/fiji/latest/fiji-linux64.zip"
     unzip fiji-linux64.zip
     # then, <path_to_fiji_or_imagej> will be "~/Fiji.app/ImageJ-linux64"
     ```

If each file in the folder is a z-stack or a RGB file already, you don't need to specify anything else. 

Otherwise, if each file is a FOV and a specific Z plane (e.g., file names like Image_00001_Z001.tif, Image_00001_Z002.tif,...), 
you can specify a regex to parse the Z plane from the file names, and a single file per FOV will be created before stitching. 
This can be specified with the argument `--join-zstack-regex`, which by default is set to an empty string
(i.e., no merging of Z planes is performed).

## (Optional) Addressing image irregularities
Most of the times, large images have uneven illumination, focus or noise. This can be challenging when doing downstream
processing on these images, like segmentation, feature extraction or quantification (i.e., on fluorescence images). 
There is a plethora of methods to address these issues (e.g., [Flatfield Correction](https://imagej.net/plugins/bigstitcher/flatfield-correction)
from *BigStitcher*, or [CARE](https://imagej.net/plugins/care), to name some). This might be highly dependent on your microscope,
imaging settings, sample type, sample width... **Always look at your images** so you can take an informed decision.

In our publication, we leveraged a [CUT](https://github.com/taesungp/contrastive-unpaired-translation) model that allowed us
to homogeneize the *style* of the whole tile-scan - that is, reduce possible biases in illumination, noise and focus - across the entire
tile-scan. You can run this by running the following command on the stitched image.

![CUT restoration showcased in the paper](../static/img/compu_image_preprocessing/preprocessing_CUT.webp){ loading=lazy }

```bash
openst image_preprocess \
     --input=<path_to_input_image> \
     --CUT \
     --CUT-model=<path_to_model> \
     --output=<path_to_output>
```

Make sure to replace the placeholders (`<...>`). For instance,
`<path_to_input_image>` is the full path and file name of the previously stitched image; `<path_to_model>`
is filename our pre-trained [CUT model](https://github.com/rajewsky-lab/openst/models/CUT.pth), and `<output_image>` 
is the path to a folder (writeable) and desired filename for the output image.

## Expected output
After running the stitching (and optionally correction algorithms), you will have a single image file per sample. This, together with
spatial transcriptomics data from the previous step, will be used in the following to align both modalities, and eventually obtain a
file that can be used for the downstream spatial, single-cell analysis.