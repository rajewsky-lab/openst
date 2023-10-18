# Preprocessing of imaging data
We'll start with the initial processing of H&E-stained tissue sections, as a prerequisite to aligning spatial transcriptomics data with tissue staining images.

## Stitching (microscope-dependent)
To begin, stitch together the tile-scan images of your H&E-stained tissue sections using the Grid/Collection stitching plugin included in [Fiji](https://imagej.net/software/fiji/downloads). This will create a composite image of the entire section.

!!! warning
     This step depends on the microscope used for imaging. In our implementation, we used a Keyence BZ-X710 inverted fluorescence phase contrast microscope, for which we provide open-source image stitching code - runs in any computer, independent of the software provided by the manufacturer. Refer to the documentation of your microscope for the stitching of tile-scans.

     If you use the same configuration described in the original open-ST publication, you can reproduce the stitching with the following code:

     ```bash
     openst image_stitch --microscope='keyence' --imagej-bin=<path_to_fiji_or_imagej> --tiles-dir=<path_to_tiles> --tiles-prefix=<to_read> --tmp-dir=<tmp_dir>
     ```

## (Optional) Style transfer
Optionally, you can perform style transfer on your tile-scan images using a custom Contrastive Unpaired Translation (CUT) model. This can help equalize the style between sections and remove artifacts.

```bash
openst image_preprocess --input=<path_to_input_image> --CUT --CUT-model=<path_to_model> --output=<path_to_output>