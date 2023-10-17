# Open-ST: Experimental Protocol

## Introduction

## Materials 

### Reagents 
??? example
        |REAGENT|SOURCE|IDENTIFIER|
        |----|----|----|
        | Dra I enzyme|NEB|Cat#R0129|
        |Alkaline Phosphatase Calf Intestinal (CIAP) enzyme|Promega|Cat#M1821|
        |Exonuclease I enzyme|NEB|Cat#M0293|
        |NaOH solution for molecular biology 10 M in H2O|Sigma|Cat#72068|
        |UltraPure™1M Tris-HCl pH7.5|Invitrogen™|Cat#15567027|
        |Tissue-Tek OCT|Sakura Finetek|Cat#4583|
        |Methanol (min. 99.8%)|Th. Geyer|Cat#1437|
        |2-propanol (min 99.9%)|Th. Geyer|Cat#1197|
        |Mayer’s Haematoxylin|Agilent Dako|Cat#S3309|
        |Bluing buffer|Agilent Dako|Cat#CS702|
        |Eosin Y, aqueous|Sigma|Cat#HT110216|
        |Pepsin from porcine gastric mucosa|Sigma|Cat#P7000|
        |20x SSC|Sigma|Cat#S6639-1L|
        |Hydrochloric Acid (HCl) 10N|AppliChem|Cat#187051|
        |BSA Molecular Biology Grade (conc. 20 mg/ml)|NEB|Cat#B9000S|
        |dNTP SET 100mM 4X1mL|Life Technologies|Cat#R0182|
        |SuperScript IV Reverse Transcriptase|Life Technologies|Cat#18090010|
        |RiboLock RNase Inhibitor|Thermo Scientific|Cat#EO0381|
        |Tris-HCl Buffer pH 8.0, 1M|Life Technologies|Cat#AM9855G|
        |Sodium chloride NaCl (5M), RNase-free|Invitrogen|Cat#AM9760G|
        |Roti®-Stock 20 % SDS ready-to-use, sterile filtered|Roth|Cat#1057.1|
        |UltraPure 0.5M EDTA, pH 8.0|Life Technologies|Cat#15575020|
        |Proteinase K (800 mU/uL)|NEB|Cat#P8107S|
        |DNA Polymerase Large Fragment exo- Klenow Fragment (3'-5' exo-)|NEB|Cat#M0212|
        |Ampure XP beads|Beckman coulter|Cat#A63881|
        |Kapa HiFi Hotstart Readymix KK2612|Roche|Cat#7958960001|
        |1.5% Agarose gel, PippinHT (300-1500 bp)|Biozym|HTC1510|
        |Qubit dsDNA HS Assay Kit|Invitrogen|Cat#Q32854|
        |High sensitivity DNA kit|Agilent|Cat#5067-4626|
        |HS RNA tapestation|Agilent|5067-5579/5580/5581|
        |Blue S'Green qPCR mix|Biozym|Cat#331416|
        |KAPA LQ Primer + Mastermix (Illumina/ LC480)|Roche|Cat#7960573001|
        |KAPA Library Quantification DNA Standards (Illumina)|Roche|Cat#7960387001|
        |NovaSeq 6000 S4 reagent kit v1.5 (35 cycles)|Illumina|Cat#20044417|

*Add oligo list*

### Equipment 

- Chemical hood (*for work with toxic chemicals, such as Trizol or methanol*)  
- Cryostat 
- Heating block (*Can also use hybridisation oven for pre-warming pepsin.*)
- Hybridization oven 
- Brightfield imaging system *add objective requirements, camera*
- PCR cycler 
- qPCR machine 
- Bluepippin or PippinHT (*alternatively, use manual agarose gel setup and DNA extraction*)
- Automated gel electrophoresis machine (eg. Tapestation or BioAnalyzer)
- Qubit fluorometer   

## 1. Capture area generation
The following section details the generation of capture areas for the open-ST protocol. By sequencing oligos, which comprise unique 32-nucleotide barcodes, appropriate adapters, and a poly-dT, we register the barcode sequences and their associated coordinates on the flow cell. 
From a NovaSeq6000 S4 flow cells ~360 capture areas sized 3x4 mm can be made.

### 1.1 Sequencing of barcoded library 

When using an Illumina® NovaSeq 6000 S4 flow cell (35 cycles), sequence the HDMI32-DraI library at a loading concentration of 200 pM. 
Sequence a single-end 37 cycle read, using Read1-DraI oligo as a custom primer. Use a custom sequencing recipe that stops the run immediately after read 1 prior to on-instrument washes. 

### 1.2 Enzymatic processing

!!! Tip 

    If bubbles occur, mark these with pen on the flow cell. Repeat reactions if many bubbles occur and ensure bubbles do not form at the same locations. 

!!! Tip 

    Use a P1000 pipette and pipette slowly to avoid the formation of bubbles. 

!!! Tip 

    For removing washes, pipette the liquid out and then blow through air using the P1000 to remove remaining liquid.

#### 1.2.1 Dra I digestion 

**DraI mix**

|Reagent|Final concentration|Volume (uL)|
|----|----|----|
|DraI|2U/uL|10|
|10X CutSmart buffer|1x|10|
|Ultrapure water||80|

1. Wash flow cell by flowing through 500 uL ultrapure water using a P1000 pipette.  
2. Add DraI mix and incubate at 37°C overnight.

#### 1.2.2 Exonucelase I digestion 

**ExoI mix**

| Reagent | Final concentration | Volume (uL) |
|----|----|----|
|ExoI|1 U/uL|5|
|10X ExoI buffer|1x|10|
|Ultrapure water||85|

1. Wash flow cell by flowing through 500 uL 80% ethanol, then ultrapure water.
2. Add Exonuclease I mix and incubate for 45 min at 37°C.
3. Wash flow cell by flowing through 500 uL ultrapure water three times.

### 1.3 Opening, denaturation and washes 

!!! Note 

    The NovaSeq6000 S4 flow cell images the top and bottom glass. Thus, keep both and take care not to break them.  

1. Remove the flow cell from its plastic encasing. 
2. Carefully score along the sides of the flow cell using a scalpell. The blade should be in one plane with the flow cell.
3. Once all sides detach, carefully seperate the two flow cell glasses. 
 
### 1.4 Breaking the flow cells into capture areas 

!!! Note
 
    Capture areas can be stored dry in the fridge for extended periods of time. We have generated libraries from prepared capture areas stored for x months. 

We have designed a cutting guide that facilitates the breaking of the flow cell into regular capture areas. 
The 3D printing file can be found here: 
*link to stl file, video*


## 2. Tissue processing and RNA quality control 

Open-ST requires the use of unfixed fresh-frozen tissue. We recommend following 10X Visium's protocol for simultaneous freezing and embedding in their [Tissue Preparation guide](https://cdn.10xgenomics.com/image/upload/v1695417744/support-documents/CG000240_Demonstrated_Protocol_VisiumSpatialProtocols_TissuePreparationGuide_RevE.pdf)
Test for RNA quality of the OCT-embedded tissue before working with a tissue. Aim for an RNA integrity number (RIN) over 6. 

## 3. Library preparation

### 3.1 Tissue sectioning

!!! Before starting 

    Place the OCT-mounted fresh frozen tissue in a cryostat for 20 minutes at the selected cutting temperature (adjusted according to tissue). Place the capture areas at room temperature. 

!!! Before starting 

    Pre-cool 100% methanol at -20°C for subsequent fixation step. 

!!! Tip 

   Trim the excess OCT surrounding the tissue to prevent folding of OCT under or over the tissue.  

!!! Warning 

    Remove the capture area from the stage, as soon as the transfer occured, to avoid re-freezing of the tissue onto the stage. 
   
1. Slice the tissue at the selected temperature at 10 μm thickness.
2. Place the capture area (room temperature) on the tissue section on the cutting stage. The tissue will melt onto the capture area.
3. Store the capture area with tissue-side up in the cryostat until fixation. 

### 3.2 Fixation 

!!! Note

    Transport samples on dry ice until placed in -20°C methanol. 

1. Fix the tissue by placing the capture area with the tissue section in a tube containing 1 mL of 100% Methanol (pre-cooled to -20°C). 
2. Incubate -20 C for 30 min.

### 3.3 H&E staining and imaging 

**Buffered Eosin (make fresh)**

| Reagent | Volume (uL) |
|----|----|
| Eosin Y | 500 |
| 0.45M Tris-Acetate buffer pH 6.0 | 500 |
|Total| 1000 |

For incubations, add enough volume to cover the capture area completely. For all washes, wash by dipping in a beaker with Ultrapure water (500 mL). 

1. Add isopropanol to the tissues on the capture area, and incubate for 1 min, then remove the solution and dry the tissues. 
2. Add hematoxylin, and incubate for 5 min 
3. Wash the FC piece with Ultrapure water until the hematoxylin dye is completely removed (10-15x). 
4. Add bluing buffer, and incubate for 2 min. 
5. Wash the FC piece with Ultrapure water, dipping 5 times. 
6. Add buffered eosin onto tissue section for 1 min.
7. Wash the slide with Ultrapure water, dipping 10-15x. 
8. Let the capture area air-dry at room temperature for 10 min. Dry the tissues completely with no residual water. 
9. Put the FC piece face-down on a coverslip (#1.5, 24x50mm) and image on brightfield with the 20x objective.

After imaging place the capture areas with tissue section face up into a multi-well gasket. 

### 3.4 Permeabilization 

We recommend performing a pilot experiment in which you compare different permeabilization conditions using a qPCR assay. We suggest comparing different pepsin incubation times (for example, 0, 15, 30, 45, and 60 min) at 37°C.    
Follow the library preparation steps until step 3.9. Earlier amplification corresponds to a higher concentration of starting material, ie. more efficient mRNA capture. If conditions amplify together, chose the shorter time or lower pepsin concentration for permeabilization of your sample.  


1. Weigh and dissolve the pepsin to have a solution with 7 U/ul pepsin in 2xSSC pH 2.5. 
2. Dilute 1:10 with 2xSSC pH 2.5 to get the final concentration of 0.7 U/ul.
3. Prewarm permeabilization mix at 37℃ several minutes before use. 
4. Incubate at 37C for X min (ex. 15 min - 30 min - 60 min) according to the tissue used.  


### 3.5 Reverse transcription 

**Reverse transcription (RT) buffer** 

|Reagent|Final concentration|Volume (ul)|
|:----|:----|:----|
|SSIV 5X rt BUFFER|1x|20|
|RNase inhibitor (40U/ul)|1U/ul|2.5|
|Ultrapure water| |77.5|
|Total| |100|

**Reverse transcription mix** 

|Reagent|Final concentration|Volume (ul)|
|:----|:----|:----|
|SSIV 5X rt BUFFER|1x|20|
|0.1M DTT|5mM|5|
|BSA (20mg/ml)|0.187mg/ml|0.93|
|10mM dNTP mix|1mM|10|
|Superscript IV (200U/ul)|6.67 U/uL|3.33|
|Ribolock(40U/uL)|1U|2.5|
|Ultrapure water| |62.33|
|Total| |100|

1. Remove the pepsin solution. 
2. Wash the capture area carefully with 100 μl RT Buffer once.
3. Add 100 uL RT mix per capture area. Incubate overnight at 42℃. 

### 3.6 Exo I digestion

**Exonuclease I mix** 

|Reagent|Final concentration|Volume (ul)|
|:----|:----|:----|
|10 x Exo I buffer|1x|10|
|Exo I|1U/ul|5|
|Ultrapure water| |85|
|Total| |100|

1. Remove the RT solution. 
2. Add 100 uL Exonuclease I mix per capture area to eliminate DNA that did not hybridize with mRNA.  
2. Incubate 45 min at 37℃. 

### 3.7 Tissue removal 

**Tissue removal mix**	

|Reagent|Final concentration|Volume (uL)|
|:----|:----|:----|
|Tris-Cl pH 8.0|100 mM|10|
|2M NaCl|200 mM|10|
|20% SDS|2%|10|
|0.1M EDTA|5 mM|5|
|Proteinase K ( 800 mU/uL)|16 mU/uL|2|
|Nuclease -free water| |63|
|Total| |100|

1. Remove the Exonuclease I mix. 
2. Add 100 μl of 1x tissue removal mix per capture area and incubate for 40 minutes at 37℃.
3. Wash as follows: 
    1. Wash the capture area with ultrapure water three times. 
    2. Wash the capture area  with 100 μl of freshly prepared 0.1N NaOH three times* (*each with 5 min incubation at room temperature). 
    3. Wash the capture area  with 100 μl of 0.1M Tris (pH7.5) three times. 
    4. Wash the capture area with 100 μl of Ultrapure water three times. 

!!! Note 

    Visually confirm that tissue removal is complete after washes have been completed. 

### 3.8 Second strand synthesis 

**Second strand synthesis mix**

|Reagent|Final concentration|Volume (uL)|
|:----|:----|:----|
|NEBuffer-2|1x|10|
|100 uM randomer|10 uM|10|
|10 mM dNTPs|1 mM|10|
|Klenow exo (-) Fragment (5 U/uL)|0.5 U/uL|10|
|Ultrapure water| |60|
|Total | |100|

1. Add 100 uL second strand synthesis mix per capture area. Incubate at 37°C for 2 h. 
2. Wash with 100 uL ultrapure water 3 times. 
3. Elute the second strand product by incubating the capture areas in 100 μl of freshly prepared 0.1 N NaOH twice for 5 min each. Recover the elutions and pool the second strand product together per sample.
4. Mix 200 μl of second strand product with 28,6 μl of 1M Tris-HCl pH 7.5. Proceed directly to next step.

Purify the 228.6 uL elution using AmpureXP beads at a ratio of 1.8 beads to 1x second strand product, following the manufacturer's instructions.
Elute the product in 82.5 uL ultrapure water.  

### 3.9 qPCR for cycle number assessment 

**qPCR mix**

|Reagent|Final concentration|Volume (ul)|
|----|----|----|
|2x Blue S#Green qPCR mix + ROX |1x|10|
|10 uM Fw primer|1 uM|2|
|10 uM Rw primer|1 uM|2|
|Second strand product| |2.5|
|ultrapure water| |3.5|
|Total| |20|

|Temperature|Time|Cycles|
|:----|:----|:----|
|95°C |3 min|1|
|95°C|30 sec|40|
|60°C|1 min| |
|72°C|1 min| |

*add example profile*


### 3.10 Library construction 

#### 3.10.1 Library amplification and purification

**Library amplification mix**

|Reagent|Final concentration|Volume (uL)|
|----|----|----|
|2x KAPA HiFi Hotstart Readymix|1x|100|
|100 uM WTA1*F primer |1 uM|2|
|100 uM WTA1*R primer|1 uM|2|
|Purified 2nd strand | |80|
|Ultrapure water| |16|
|Total| |200|


|Temperature|Time|Cycles|
|:----|:----|:----|
|95°C |3 min|1|
|95°C|30 sec|To be determined|
|60°C|1 min| |
|72°C|1 min| |
|72°C|2 min|1|
|4°C|hold|

1. Prepare the library amplification mix per sample. 
2. Split wach sample mix into four PCR tubes, each with 50 uL volume. 
3. Run the PCR with the cycle number determined previously (3.9). 
4. Pool the 200 uL PCR product per sample and purify using AmpureXP beads at a 1:1 ratio of beads PCR product, following the manufacturer's instructions.
5. Elute the PCR product in 20 uL ultrapure water. 

#### 3.10.2 Size selection 

Perform size selection of your sample to obtain fragments 350 - 1100 bp. Use the Bluepippin or PippinHT 1.5% agarose gel and follow the manufacturer's instructions. 
Measure the concentration of your size-selected product using the Qubit dsDNA quanitification kit. 

## 4. Sequencing



