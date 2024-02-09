# Preprocessing of imaging data

For this dataset, we archived the stitched tile-scan image (a link will be available soon). 
This single image was generated from multiple, independently imaged tiles, by leveraging `openst image_stitch`. So, 
you don't need to use this command, since we already provide the stitched image. Anyway, let us know
if you want access to this tile images, in case you want to try.

As well, the imaging from this dataset did not require any further postprocessing prior to segmentation,
as visual inspection of the images did not reveal any strong illumination or focus biases. 

*(provide a screenshot of the image here)*

Therefore, feel free to download the data from the link above, and continue to the next section.
Then, browse back to the `openst_e13_demo` folder, create a folder `alignment`, and save the
`e13_brain_staining.tif` file inside the `alignment` folder. You will have the following file structure:

```sh
/home/user
|-- openst_e13_demo
|   |-- data
|   |   `-- fastq
|   |-- spacemake
|   |   `-- ...
|   |-- bins
|   `-- alignment
|       `-- e13_brain_staining.tif
```