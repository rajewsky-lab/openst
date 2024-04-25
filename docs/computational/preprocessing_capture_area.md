# Preprocessing capture area library
After [sequencing the Open-ST capture areas](../experimental/capture_area_generation.md), you will get basecall files in `bcl` format, or raw reads in `fastq`
format (see [sequence file formats](https://www.illumina.com/informatics/sequencing-data-analysis/sequence-file-formats.html)
from Illumina's website).

We have designed a simple computational workflow to transform the raw `bcl` or `fastq` files
from the sequencing of the barcoded library into table-like files (`csv`, or `tsv`) with contain the following
information:

|cell_bc|x_pos|y_pos|
|----|----|----|
|CGCGAGGGGAAAATGGGGACTAGCG|6343|1016|
|GGTCCCGTCCAAGAAGTAAATCGAA|9272|1016|
|...|...|...|

Where `cell_bc` is the 32 nucleotide-long spatial barcode, and `x_pos`/`y_pos` are 2D spatial coordinates
of a specific **tile** in the capture area (see below). 

Before diving into the code, let's clarify some of the
[terms](#flow-cell-related-terms) that are specific to using Illumina flow cells as capture areas. We quote from 
[Illumina's documentation](https://support-docs.illumina.com/IN/NextSeq_550-500/Content/IN/NextSeq/FlowCell_Tiles_fNS.htm)

### Flow cell-related terms
!!! info "Tiles"
    "Small imaging areas on the flow cell defined as the field of view by the 
    camera. The total number of tiles depends on the number of lanes, swaths, and surfaces that are imaged on 
    the flow cell, and how the cameras work together to collect the images."

!!! info "Lane"
    "A physical channel with dedicated input and output ports."

!!! info "Top/bottom"
    "The flow cell is imaged on two surfaces, the top and bottom.
    The top surface of 1 tile is imaged, then the bottom surface of the same tile is imaged before moving to the next tile."

!!! info "Swath"
    "A column of tiles in a lane."

### Retrieve spatial barcodes coordinates for one tile
The `x_pos` and `y_pos` coordinates from the table above are given for each tile, separately. This information is
encoded in the `bcl` and `fastq` files. To obtain per-tile barcodes and coordinates, run the following code: 


```sh
openst barcode_preprocessing \
    --fastq-in <fastq_of_tile> \
    --tilecoords-out <out_path> \
    --out-suffix <out_suffix> \
    --out-prefix <out_prefix> \
    --crop-seq <len_int> \
    --rev-comp \
    --single-tile
```

Make sure to replace the placeholders:
`<fastq_of_tile>` to the `fastq` file of a specific tile; `<out_path>` where the table-like
files will be written; `<out_suffix>` and `<out_prefix>` are suffixes and prefixes that are added to the tile file names;
`<len_int>` from the `--crop-seq` argument is a string in the [Python slice](https://docs.python.org/3/tutorial/datastructures.html) format
(e.g., 2:32 will take nucleotides 2nd until 32th of the sequence in the `fastq` file); `--rev-comp` is provided whether the barcode sequences
must be written into the `csv` as their reverse-complementary; `--single-tile` argument is provided when the `fastq` file only contains data for
a single tile (**our recommendation**).

### Retrieve spatial barcodes coordinates for all tiles
Above you generated a single tile coordinate. To process all tiles from a flow cell (in parallel), you can run the 
following snippets for Linux, assuming you have access to the basecalls folder.

First create a `lanes_and_tiles.txt` file:

```sh
cat RunInfo.xml | grep "<Tile>" | sed 's/ *<Tile>//' | sed 's/<\/Tile>//' | sed 's/^[ \t]*//;s/[ \t]*$//' > lanes_and_tiles.txt
```

where `RunInfo.xml` is a file contained in the basecalls directory.

Then, run demultiplexing and conversion to `fastq`
simultaneously to generating the barcode spatial coordinate file:

```sh
cat lanes_and_tiles.txt | xargs -n 1 -P <parallel_processes> -I {} \
    sh -c 'bcl2fastq -R <bcl_in> --no-lane-splitting \
                -o <bcl_out>/"{}" --tiles s_"{}"; \

            openst barcode_preprocessing \
                --fastq-in <bcl_out>/{}/Undetermined_S0_R1_001.fastq.gz \
                --tilecoords-out <out_path> \
                --out-suffix .txt \
                --out-prefix <out_prefix>"{}" \
                --crop-seq <len_int> \
                --rev-comp \
                --single-tile'
```

Again, make sure to replace the placeholders: `<bcl_in>` and `<bcl_out>` are the directories where the
basecall files are contained and where the converted output `fastq` files will be saved; The rest of
arguments have the same meaning as above. If you generated full `fastq` yourself, you can adapt the command
above to remove the call to `bcl2fastq`.

## Expected output

After running all the steps of this section, you will have a folder with many `*.txt.gz` files
containing the spatial coordinates of flow cell tiles (you will only need to generate this once per flow cell),
in the tab format described above.