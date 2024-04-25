# Sequencing of Open-ST library

## Quantification
We recommend quantifying your libraries for sequencing using the KAPA Library Quantification Kit.

## Loading and sequencing
The optimal loading concentration depends on the sequencer used. For the Illumina® **NovaSeq 6000** we obtained optimal clustering at a loading concentration of **130 pM**. For the Illumina® **NextSeq 2000** sequencing system we recommend a loading concentration of **650 pM**.   

Moreover, Illumina [suggests a minimum spike-in of 1% PhiX](https://support.illumina.com/content/dam/illumina-support/documents/documentation/system_documentation/novaseq/novaseq-6000-denature-dilute-libraries-guide-1000000106351-03.pdf) into the pool as a quality control for cluster generation, sequencing, and alignment. 

In our setup, the following read lengths were used:

|Read|Cycles|
|----|----|
|Read 1|28-32|
|Index 1|8|
|Index 2|NA|
|Read2|90+|

## Expected (data) output
Either when using your own sequencing equipment or relying on a sequencing facility, you will get access
to (most likely) already [demultiplexed](https://knowledge.illumina.com/software/general/software-general-troubleshooting-list/000005982)
`fastq` files; otherwise, you can get access to the *raw* basecall files in `bcl` [format](https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/ToolsiBCL_fDG.htm).

Either of these files shall be used as the input for spacemake [later](../computational/preprocessing_openst_library.md#processing-of-the-Open-ST-library).
