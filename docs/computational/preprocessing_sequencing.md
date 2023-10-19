# Preprocessing sequencing data
After sequencing, you will get basecall files in `bcl` format, or raw reads in `fastq`
format (see [sequence file formats](https://www.illumina.com/informatics/sequencing-data-analysis/sequence-file-formats.html)
from Illumina's website). In the open-ST experimental protocol, there are two points at which
sequencing is performed: (1) for obtaining the [sequences of the barcoded library](experimental.md#1-1-sequencing-of-barcoded-library),
and (2) for obtaining the [transcriptomic sequences of open-ST library](experimental.md#4-sequencing-of-open-st-library).

We have developed computational tools for processing each of these sequencing rounds.
By processing (1), you will get a *database* of barcodes and their spatial locations.
This will used during the automated processing of (2), such that the transcriptome is mapped back to space.

## Processing the sequencing of barcoded library
We have designed a simple computational workflow that allows to transform the `bcl` or `fastq` files
from the sequencing of the barcoded library into table-like files (`csv`, or `tsv`) that contain the following
structure:

|cell_bc|x_pos|y_pos|
|----|----|----|
|CGCGAGGGGAAAATGGGGACTAGCG|6343|1016|
|GGTCCCGTCCAAGAAGTAAATCGAA|9272|1016|
|...|...|...|

Where `cell_bc` is the 32 nucleotide-long spatial barcode, and `x_pos`/`y_pos` are the 2d spatial coordinates
of a specific **tile** in the capture area (see below). Before diving into the code, let's clarify some of the
[terms](#flow-cell-related-terms) that are specific to using Illumina flow cells as capture areas. We quote from 
[Illumina's documentation](https://support-docs.illumina.com/IN/NextSeq_550-500/Content/IN/NextSeq/FlowCell_Tiles_fNS.htm)

### Flow cell-related terms
!!! info
    **Tiles**: "small imaging areas on the flow cell defined as the field of view by the 
    camera. The total number of tiles depends on the number of lanes, swaths, and surfaces that are imaged on 
    the flow cell, and how the cameras work together to collect the images."

!!! info
    **Lane**: "a physical channel with dedicated input and output ports."

!!! info
    **Top/bottom**: "the flow cell is imaged on two surfaces, the top and bottom.
    The top surface of 1 tile is imaged, then the bottom surface of the same tile is imaged before moving to the next tile."

!!! info
    **Swath**: "a column of tiles in a lane."

### Computing barcodes and spatial coordinates of all tiles
The `x_pos` and `y_pos` coordinates from the table above are given for each tile, separately. This information is
encoded in the `bcl` and `fastq` files. To obtain per-tile barcodes and coordinates, run the following code: 


```sh
openst barcode_preprocessing \
    --in-fastq <fastq_of_tile> \
    --out-path <out_path> \
    --out-suffix <out_suffix> \
    --out-prefix <out_prefix> \
    --crop-seq <len_int> \
    --rev-comp \
    --single-tile
```

Make sure to replace the placeholders (`<...>`). For instance,
`<fastq_of_tile>` is the full path to the `fastq` file of a specific tile; `<out_path>` is the folder where the table-like
files will be written; `<out_suffix>` and `<out_prefix>` are suffixes and prefixes that are added to the tile file names;
`<len_int>` from the `--crop-seq` argument is a string in the [Python slice](https://docs.python.org/3/tutorial/datastructures.html) format
(e.g., 2:32 will take nucleotides 2nd until 32th of the sequence in the `fastq` file); `--rev-comp` is provided whether the barcode sequences
must be written into the `csv` as their reverse-complementary; `--single-tile` argument is provided when the `fastq` file only contains data for
a single tile (**our recommendation**).

The code above will generate a file in `<out_path` per tile. Only a single fastq file can be provided at a time via `--in-fastq`. To
process this in parallel, you can run the following snippets (in Linux, assuming you start from the `fastq` files). We assume that
you have a file `lanes_and_tiles.txt`, that contains the tile identifiers that you want to process; you can generate this file with:

```sh
cat RunInfo.xml | grep "<Tile>" | sed 's/ *<Tile>//' | sed 's/<\/Tile>//' | sed 's/^[ \t]*//;s/[ \t]*$//' > lanes_and_tiles.txt
```

where `RunInfo.xml` is a file contained in the basecalls directory. *We don't endorse parsing xml like this, but
this code snippet works* 🙈. Then, you can process various `fastq` files in the basecalls directory as follows:

```sh
cat lanes_and_tiles.txt | xargs xargs -n 1 -P <parallel_processes> -I {} \
    sh -c 'openst barcode_preprocessing \
                --in-fastq <fastq_dir>/{}/Undetermined_S0_R1_001.fastq.gz \
                --in-fastq <fastq_of_tile> \
                --out-path <out_path> \
                --out-suffix .txt \
                --out-prefix <out_prefix>"{}" \
                --crop-seq <len_int> \
                --rev-comp \
                --single-tile'
```

Make sure to replace the placeholders (`<...>`). For instance,
`<parallel_processes>` is the number of parallel processes that will be spawned (recommended: less than the number
of cores in your machine); `<fastq_dir>` is the subdirectory of the basecalls directory where `fastq` files are contained; 
`<out_prefix>` is the prefix added to the file names (e.g., fc_1 as an internal unique identifier for a flow cell,
so you can keep track when having more than one flow cell).

Otherwise, if you start from `bcl` files (raw basecalls), you can run demultiplexing and conversion to `fastq`
simultaneously to generating the barcode spatial coordinate file:

```sh
cat lanes_and_tiles.txt | xargs xargs -n 1 -P <parallel_processes> -I {} \
    sh -c 'bcl2fastq -R <bcl_in> --no-lane-splitting \
                -o <bcl_out>/"{}" --tiles s_"{}"; \

            openst barcode_preprocessing \
                --in-fastq <bcl_out>/{}/Undetermined_S0_R1_001.fastq.gz \
                --in-fastq <fastq_of_tile> \
                --out-path <out_path> \
                --out-suffix .txt \
                --out-prefix <out_prefix>"{}" \
                --crop-seq <len_int> \
                --rev-comp \
                --single-tile'
```

Again, make sure to replace the placeholders (`<...>`). Now, `<bcl_in>` and `<bcl_out>` are the directories where the
basecall files are contained and where the converted output `fastq` files will be saved; The rest of
arguments have the same meaning as above.

## Processing of the open-ST library
The transformation of raw sequencing data into spatially-mapped expression matrices was carried out utilizing 
[spacemake](https://spacemake.readthedocs.io/en/latest/) (see also in [GitHub](https://github.com/rajewsky-lab/spacemake)), 
an automated pipeline designed for the preprocessing, alignment, and quantification of single-cell and spatial transcriptomics data.

We refer to the [official documentation](https://spacemake.readthedocs.io/en/latest/) for a complete tutorial on how to install and
run spacemake. In summary, the user needs to specify a configuration file (by default, we included run modes that are compatible
with open-ST), and a project configuration, containing the locations to the `fastq` files, to the **spatial barcode coordinate files**
(generated in the [previous step](#processing-the-sequencing-of-barcoded-library)),
and other metadata. Then, spacemake automatically processes these data into [h5ad](https://anndata.readthedocs.io/en/latest/fileformat-prose.html)
files, which basically contain a matrix of barcodes (rows) and genes (columns), with associated metadata and spatial coordinates. 

As well, several HTML reports are generated by spacemake, such that the user can quickly and visually
assess the quality of samples prior to any downstream processing and analysis.

## Global spatial coordinates: *tile* stitching
The `x_pos` and `y_pos` coordinates from the [table above](#processing-the-sequencing-of-barcoded-library) are given relative to each individual tile, not absolute to the 
whole extension of the flow cell - which consists of many *tiles* arranged in *lanes* and *swaths*, see above. 
This manner, we provide code in `openst` and `spacemake` that allows to compute the global coordinates in the 
case of having a sample whose placement has cover more than one tile
during the open-ST library preparation. 

When specifying a `run_mode` without meshing and with a `puck_collection`
(see [spacemake documentation](https://spacemake.readthedocs.io/en/latest/config.html),
under the section *Configure run_modes*), spacemake will automatically generate a single file out of 
all individual (per-tile) `h5ad` files, with the suffix `_puck_collection`. This file must be used for all
subsequent steps of this tutorial, not the individual `h5ad` files.

Otherwise, if these `puck_collection` files were not automatically generated, we provide code within the
`openst` package with that same functionality. This must be ran one sample at a time, and requires two types
of input file: (1) all the `h5ad` files generated for a single sample (can be specified implicitly, via a wildcard in Linux,
or explicitly as a list separated by spaces); and (2) a coordinate system file, that specifies that is the relative offset,
in units equivalent to `x_pos` and `y_pos`, of the tiles respect to their *column* and *swath*. These coordinate system files
have a standardized format, and are provided in the [spacemake] and [openst] repos for the Illumina® NovaSeq 6000 S4 flow cell.
If you have used other flow cell, you might need to generate these files following the same convention (feel free to open an
issue in our [repo], and we can generate it for you).

To manually create 'puck_collection' files, you can run the following in a terminal:

```sh
openst spatial_stitch \
    --tiles <space_separated_list_or_wildcards_to_h5ad> \
    --tile-coordinates <path_to_coordinate_system> \
    --output <output_puck_collection_h5ad>
```

This program has additional arguments that are explained when running `openst spatial_stitch --help`. Make sure to replace
the placeholders (`<...>`); for example, `<space_separated_list_or_wildcards_to_h5ad>` is a space-separated list or a
implicit (wildcards) path to all h5ad of tiles for a single sample, from spacemake (at the automatically generated `dge` folder).
The `<path_to_coordinate_system>` is a path to the `csv` file containing the relative offsets of tiles; and, `<output_puck_collection_h5ad>`
is the name (full or relative path) of the file that will be generated.

## Expected output

After running all the steps of this section, you will end up with two types of file: (1) the spatial coordinates of flow cell
tiles (you will only need to generate this once per flow cell); and (2) one `h5ad` file per sample, containing the gene expression and
spatial coordinates of each barcoded spot. In the following sections, we will use the file (2) and the 
[previously preprocessed images](computational/preprocessing_imaging.md)