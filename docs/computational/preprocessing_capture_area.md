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

Where `cell_bc` is the 25 nucleotide-long spatial barcode, and `x_pos`/`y_pos` are 2D spatial coordinates
of a specific **tile** in the capture area (see below). 

Before diving into the code, let's clarify some of the
[terms](#flow-cell-related-terms) that are specific to using Illumina flow cells as capture areas. We quote from 
[Illumina's documentation](https://support-docs.illumina.com/IN/NextSeq_550-500/Content/IN/NextSeq/FlowCell_Tiles_fNS.htm)

## Flow cell-related terms
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

## Generating barcode-to-coordinate map

For each flow cell, we generate plain text files with three columns: `cell_bc`, `x_pos`, and `y_pos`.
These files are later used by `spacemake` to reconstruct the spatial coordinates from transcriptomic libraries. 
This process is performed only once per barcoded flow cell.

To create the barcode-to-coordinate map for all tiles, use the following command:

```sh
openst flowcell_map \
    --bcl-in /path/to/fc/bcl \
    --tiles-out /path/to/fc_tiles \
    --crop-seq 5:30 \  # for default Open-ST sequencing recipe
    --rev-comp
```

Make sure to specify the arguments:
- `--crop-seq`: Use a compatible Python slice (e.g., 5:30 will take 25 nucleotides, from the 6th to the 30th from the input reads)
- `--rev-comp`: After cropping the sequences, will compute and store the reverse complement of the barcode sequences

This command will create as many barcode-to-coordinate compressed text files as there are tiles in the flow cell under the folder `/path/to/fc_tiles`

## (Optional) Retrieve spatial barcodes coordinates for one tile
The `x_pos` and `y_pos` coordinates from the table above are given for each tile, separately. This information is
encoded in the `bcl` and `fastq` files. To obtain per-tile barcodes and coordinates, run the following code: 


```sh
openst barcode_preprocessing \
    --fastq-in /path/to/tile.fastq \
    --tilecoords-out /path/to/fc_tiles \
    --out-prefix fc_1_ \
    --crop-seq 5:30 \
    --rev-comp \
    --single-tile
```

Make sure to replace the placeholders.
`/path/to/tile.fastq` to the `fastq` file of a specific tile; `/path/to/fc_tiles` where the table-like
files will be written; `--out-prefix` (and `--out-suffix`) are prefixes and suffixes that are added to the tile file names;
`--crop-seq 5:30` is a [Python slice](https://docs.python.org/3/tutorial/datastructures.html)
(e.g., 5:30 will take nucleotides 6th until 30th of the sequence in the `fastq` file); `--rev-comp` is provided whether the barcode sequences
must be written into the `csv` as their reverse-complementary, **after cropping**; `--single-tile` argument is provided when the `fastq` file only contains data for
a single tile (**our recommendation**).

## Expected output

After running all the steps of this section, you will have a folder `/path/to/fc_tiles` with `*.txt.gz` files
containing the spatial coordinates of flow cell tiles. **You only need to generate this once per flow cell.**