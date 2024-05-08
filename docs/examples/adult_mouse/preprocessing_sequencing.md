# Preprocessing transcriptomic library

After sequencing, we proceed with the preprocessing of the data, to go from raw reads to
transcriptomic information mapped to the mouse genome, in space.

## Demultiplexing
We got basecall files in `bcl` format from our sequencing facility.

We used `bcl2fastq` for demultiplexing, using this [sample sheet](../../static/examples/adult_mouse/sample_sheet.csv).
We used the conda environment where we installed `spacemake` (see [instructions
on how to install spacemake](https://spacemake.readthedocs.io/en/latest/install.html)), and ran the following commands:

```sh
(base) user@computer:~$ bcl2fastq \
    --no-lane-splitting \
    --runfolder-dir /openst/data/0_basecalls/230616_VH01346_22_AAC2LVVHV \
    -o /openst/data/1_spacemake_mouse/demultiplexed_data \
    --sample-sheet /openst/data/0_sample_sheets/230616_NR_FC_ST_72_76_AT_01.csv
```

We obtained `fastq` files that will be used for the rest of the pipeline, for [Read 1](http://bimsbstatic.mdc-berlin.de/rajewsky/openst-public-data/adult_mouse_hippocampus_R1_001.fastq.gz) and [Read 2](http://bimsbstatic.mdc-berlin.de/rajewsky/openst-public-data/adult_mouse_hippocampus_R2_001.fastq.gz).
Once you download these files, you can move them anywhere in your filesystem. We assume that you have opened a terminal,
and you have browsed to your home directory. From there, create a folder `openst_adult_demo`; browse inside, and create
another folder `data`. Then, copy the folder with the `fastq` files in here. You should have a structure like:

```sh
/home/user
|-- openst_adult_demo
|   `-- data
|       `-- fastq
```

## Transcriptomic & spatial mapping with `spacemake`
First of all, intialize the conda environment for `spacemake`
```sh
(base) user@computer:~$ conda activate spacemake
(spacemake) user@computer:~$
```

### Initialize
Create two folders inside your `openst_adult_demo` folder, called `spacemake` and `bins`, so you will have:

```sh
/home/user
|-- openst_adult_demo
|   |-- data
|   |   `-- fastq
|   |-- spacemake
|   `-- bins
```

Download the [DropSeq tools](https://github.com/broadinstitute/Drop-seq/releases/download/v2.5.4/Drop-seq_tools-2.5.4.zip),
decompress it, and put it inside the `bins` subdirectory.

Then, following the [spacemake *Quick start guide*](https://spacemake.readthedocs.io/en/latest/quick-start/index.html#open-st-quick-start),
browse to the `spacemake` directory you just created in the `openst_adult_demo` folder, and run the initialization:

```sh
(spacemake) user@computer:~$ cd /home/user/openst_adult_demo/spacemake
(spacemake) user@computer:/home/user/openst_adult_demo/spacemake$ spacemake init
    --dropseq_tools /home/user/bins/Drop-seq_tools-2.5.1
```

### Configure

As `spacemake` comes with no default value for species, before anything can be done, a new species has to be added.
In this case, we add mouse; you will need to download the correct `fa` and `gtf` files. For instance, you can download the
mouse genome from [gencode](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M30/GRCm39.genome.fa.gz),
as well as the [annotation](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M30/gencode.vM30.annotation.gtf.gz).

Then, you need to run the following commands (remember, in the same `spacemake` folder as before, with the `spacemake` conda environment;
we are going to omit `(spacemake) user@computer:/home/user/openst_adult_demo/spacemake$` for brevity).

```sh
spacemake config add_species \
   --name mouse \
   --reference genome \
   --sequence GRCm39vM30.genome.fa \
   --annotation gencodevM30.annotation.gtf

spacemake config add_species \
   --name mouse \
   --reference rRNA \
   --sequence mouse.rRNA.fa

spacemake config add_species \
   --name mouse \
   --reference phiX \
   --sequence phiX.fa
```

!!! note
    The `.fa` and `.gtf` files for mouse are available for http download under the [example datasets](../datasets.md) page.
    For instance, you can run:

    ```sh
    # for the mouse genome sequence
    wget "http://bimsbstatic.mdc-berlin.de/rajewsky/openst-public-data/genomes/GRCm39vM30.genome.fa"
    # etc...
    ```

### Add sample

Now you need to add the sample data and metadata to `spacemake`. For this, you will also need the puck (tile) barcode files, which [can be
generated](../../computational/preprocessing_capture_area.md#computing-barcodes-and-spatial-coordinates-of-all-tiles) with the `openst` package.

For simplicity, we provide the [tile barcode files](https://bimsbstatic.mdc-berlin.de/rajewsky/openst-public-data/adult_hippocampus_tiles.tar.xz) that are related to this sample, as well as the [coordinate system 
for the Illumina flow cell](https://bimsbstatic.mdc-berlin.de/rajewsky/openst-public-data/fc_2_coordinate_system.csv) that was used to generate the capture area of this experiment.

When downloading the tile barcode files, create a folder under `openst_adult_demo/data` called `tiles`. Move the files of the tile barcode files
into this folder. Also, move the coordinate file to the `puck_data` folder in the `spacemake` directory.

Remember! You need to be in the `/home/user/openst_adult_demo/spacemake` directory (or similar, depending on what you created);
then run the following commands:

```sh
# downloading R1 and R2 sequences
wget "http://bimsbstatic.mdc-berlin.de/rajewsky/openst-public-data/adult_mouse_hippocampus_R1_001.fastq.gz"
wget "http://bimsbstatic.mdc-berlin.de/rajewsky/openst-public-data/adult_mouse_hippocampus_R2_001.fastq.gz"

# downloading the spatial barcode sequences
wget "http://bimsbstatic.mdc-berlin.de/rajewsky/openst-public-data/adult_hippocampus_tiles.tar.xz"
tar -xvf adult_hippocampus_tiles.tar.xz

spacemake projects add_sample \
    --project_id openst_demo \
    --sample_id openst_demo_adult_mouse \
    --R1 adult_mouse_hippocampus_R1_001.fastq.gz \
    --R2 adult_mouse_hippocampus_R2_001.fastq.gz \
    --species mouse \
    --puck openst \
    --run_mode openst \
    --barcode_flavor openst \
    --puck_barcode_file adult_hippocampus_tiles/*.txt.gz \
    --map_strategy "bowtie2:phiX->bowtie2:rRNA->STAR:genome:final"
```

You can specify the coordinate system by modifying the `openst` run mode in the `config.yaml` file that is created
after you run the `spacemake init` command (see above). Modify the following lines from this:

```yaml
openst:
    coordinate_system: puck_data/openst_coordinates.csv
    spot_diameter_um: 0.6
    width_um: 1200
```

into this:

```yaml
openst:
    coordinate_system: puck_data/fc_2_coordinate_system.csv
    spot_diameter_um: 0.6
    width_um: 1200
```

You can download the coordinate system file from the Open-ST website, for example:

```sh
# download the coordinate system
wget "http://bimsbstatic.mdc-berlin.de/rajewsky/openst-public-data/fc_2_coordinate_system.csv"
cp fc_2_coordinate_system.csv puck_data/.
```

### Run
That's all you need to configure! Now, you can run spacemake with the following:

```sh
spacemake run --cores 32
```

You can modify the number of `--cores` depending on your local machine; also, you can specify additional
arguments to `spacemake run` - refer to the official documentation.

## Expected output

Once `spacemake` finishes, you will see that several folders and files have been created under `projects`
(inside the `spacemake` directory). For example, you can check the QC reports in your web browser by opening the
file at `projects/openst_demo/processed_data/openst_demo_adult_mouse/illumina/complete_data/qc_sheets/qc_sheet_openst_demo_adult_mouse_fc_1_puck_collection.html`,
giving you a first impression of what's the quality of spatial mapping, amount of transcripts and genes per barcoded spot, and others.

Importantly, you will find files in the directory `projects/openst_demo/processed_data/openst_demo_adult_mouse/illumina/complete_data/dge`
with the name `dge.all.polyA_adapter_trimmed.mm_included.spatial_beads_*.h5ad` (where `*` is a wildcard). These files, and not the ones containing 
the words `mesh`, `hexagon` or `circle` are the ones that will be used later to perform the pairwise alignment with imaging data, and to
later reconstruct the cell-by-gene matrix.

If you specified options for *meshing* in the `run_mode`, there will be a file containing keywords `puck_collection` and `mesh`, `hexagon` or `circle`.
This contains *approximate* cell-by-gene information, as the transcripts are aggregated by a regular lattice and not by the true spatial arrangement of
cells. This might be already enough for some analyses. 

Anyway... keep going with the tutorial if you want to unleash the full potential of Open-ST.
