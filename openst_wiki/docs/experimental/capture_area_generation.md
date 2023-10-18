# Capture area generation
The following section details the generation of capture areas for the open-ST protocol. 

By sequencing oligos, which comprise unique 32-nucleotide barcodes, appropriate adapters,
and a poly-dT, we register the barcode sequences and their associated coordinates on the flow cell.

For instance, you can get ~360 capture areas sized 3x4 mm from a single **Illumina® NovaSeq 6000 S4** flow cell.

## Sequencing of barcoded library 

When using an **Illumina® NovaSeq 6000 S4** flow cell (35 cycles), sequence the HDMI32-DraI library
(see in [Oligonucleotides](getting_started.md)) at a loading concentration of 200 pM. 

Sequence a single-end 37 cycle read, using Read1-DraI oligo as a custom primer.
Use a custom sequencing recipe that stops the run immediately after read 1 prior to on-instrument washes. 

### Expected (data) output
Either when using your own sequencing equipment or relying on a sequencing facility, you will get access
to (most likely) already [demultiplexed](https://knowledge.illumina.com/software/general/software-general-troubleshooting-list/000005982)
`fastq` files; otherwise, you can get access to *raw* basecall files in `bcl` [format](https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/ToolsiBCL_fDG.htm).

Either of these files shall be used as the input for `openst` [later](../computational/preprocessing_sequencing.md#processing-the-sequencing-of-barcoded-library),
to create a database of barcode sequences and their spatial locations.

## Enzymatic processing

!!! Tip "Tips" 

    - If bubbles occur, mark these with pen on the flow cell. Repeat reactions if many bubbles occur and ensure bubbles do not form at the same locations. 
    - Use a P1000 pipette and pipette slowly to avoid the formation of bubbles. 
    - For removing washes, pipette the liquid out and then blow through air using the P1000 to remove remaining liquid.

### Dra I digestion 

**DraI mix**

|Reagent|Final concentration|Volume (uL)|
|----|----|----|
|DraI|2U/uL|10|
|10X CutSmart buffer|1x|10|
|Ultrapure water||80|

1. Wash flow cell by flowing through 500 uL ultrapure water using a P1000 pipette.  
2. Add DraI mix and incubate at 37°C overnight.

### Exonucelase I digestion 

**ExoI mix**

| Reagent | Final concentration | Volume (uL) |
|----|----|----|
|ExoI|1 U/uL|5|
|10X ExoI buffer|1x|10|
|Ultrapure water||85|

1. Wash flow cell by flowing through 500 uL 80% ethanol, then ultrapure water.
2. Add Exonuclease I mix and incubate for 45 min at 37°C.
3. Wash flow cell by flowing through 500 uL ultrapure water three times.

## Opening, denaturation and washes 

!!! Note 

    The NovaSeq6000 S4 flow cell images the top and bottom glass. Thus, keep both and take care not to break them.  

1. Remove the flow cell from its plastic encasing. 
2. Carefully score along the sides of the flow cell using a scalpell. The blade should be in one plane with the flow cell.
3. Once all sides detach, carefully seperate the two flow cell glasses. 
 
## Breaking the flow cells into capture areas 

!!! Note
 
    Capture areas can be stored dry in the fridge for extended periods of time. We have generated libraries from prepared capture areas stored for x months. 

We have designed a cutting guide that facilitates the breaking of the flow cell into regular capture areas. 
We provide the 3D model of the cutting guide as a printable [stl](../open-ST_cutting_guide_3D_design.stl) file.
If you don't have a 3D printer, you can check for [3D printing services near you](https://www.google.com/search?q=3d+printing+service+near+me) -
they will help you in this endeavor 🤗.


Once you have the tool, refer to the video below where we explain how to use it:

---

:fontawesome-brands-youtube:{ style="color: #EE0F0F" }
__[open-ST: breaking the flow cell into capture areas]__ by @danilexn – :octicons-clock-24:
10m – Learn how to break an Illumina® NovaSeq 6000 S4 flow cell into capture areas 
using our 3D-printable cutting guide.

  [open-ST: breaking the flow cell into capture areas]: https://www.youtube.com

---