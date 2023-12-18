# FAQs 

Do you have questions about the experimental or computational aspects of open-ST?
We do our best to answer all of your questions on this page. If you can't find your question 
below, ask it on our [discussion board]!

  [discussion board]: https://github.com/rajewsky-lab/openst/discussions/new/chooses


## Tissue handling and sectioning

[__What is the recommended freezing and embedding process?__](#recommended-freezing-embedding){ #recommended-freezing-embedding }

Open-ST requires the use of unfixed fresh-frozen tissue, embedded in OCT. An optimal embedding process is required to avoid the formation of freezing artifacts. The best results are achieved when tissue is frozen in a fast and uniform way.
We recommend following 10X Visium's protocol for simultaneous freezing and embedding in their [Tissue Preparation guide](https://cdn.10xgenomics.com/image/upload/v1695417744/support-documents/CG000240_Demonstrated_Protocol_VisiumSpatialProtocols_TissuePreparationGuide_RevE.pdf)

However, is also possible to embed already snap-frozen tissue in OCT using powdered dry ice, or an isopentane bath in liquid nitrogen or dry ice. In this case, be aware of:

- Once tissue is frozen avoid melting of the tissue when embedding in OCT.
- Do not immerge the tissue directly into liquid nitrogen (it will burn the edges!).
- Keep the embedding mold bases fully covered by dry ice or the cold isopentane bath to allow homogeneous freezing.


[__What tissue section thickness is recommended?__](#tissue-thickness){ #tissue-thickness }

We section the tissue at a thickness of 10 um. This can be adapted if neccessary, but permeabilization time may change consequently. Thicker sectioning may also increase contamination from the cytoplasm of other cells in the z-plane.

[__What cutting temperature is recommended?__](#cutting-temperature){ #cutting-temperature}

The optimal cutting temperature depends on the tissue being sectioned. As a starting point, page 91 in this [reference](../387928-EP_05 - NX70 Op Guide.pdf) by Epredia is helpful.

## Sample handling

[__Are technical replicates required?__](#technical-replicates){ #technical-replicates }

As we have observed a high reproducibility across technical replicates (consecutive sections), technical replicates are not essential.

[__How many sections (conditions/sample types) can be processed at once?__](#samples-process-at-once){ #samples-process-at-once }
	
Multiple samples can be processed at once with the use of a multiwell chamber hybridization cassette
(such as this one by [ArrayJet](https://arrayjet.com/consumables/proplate-multi-well-chambers-16-well-slide-module) (product code: 206862)) (ArrayJet, accessed 13.11.2023). 

This cassette has sixteen 7 x 7 mm wells, each fitting one capture area; thus, allowing 16 samples or conditions to be processed simultaneously.
We recommend handling a maximum of 15 capture areas per person for protocol steps 3.1 to 3.5 (until overnight incubation for reverse transcription).

Moreover, open-ST libraries are indexed on the p7-adapter side, allowing multiplexing of samples within one NGS run. 

[__What are the best practices to avoid cross-contamination?__](#best-practices-x-contamination){ #best-practices-x-contamination }

In order to avoid any RNA contamination is important to wipe down the cryostat and any tools used (brushes, tweezers) with 80%-100% ethanol before sectioning. Additionaly, change the blade and wipe down the cryostat and tools in between sectioning different samples. 

[__What pepsin incubation times should I test?__](#pepsin-incubation){ #pepsin-incubation }

It is important to set the permeabilization condition for each tissue type. 
We recommend to test at least a range including 15 min, 30 min, and 60 min with two different concentration (0.7 and 1.4 U/uL). 

It is preferable to chose the minimum incubation time/ enzyme concentration that
gives the maximum RNA capture (see [Permeabilization](experimental/library_preparation.md#permeabilization)).

[__What tissues have already been tested with open-ST?__](#tested-tissues){ #tested-tissues }

Several tissues have been tested using open-ST, including *in-vitro* 3D-cultures.
Here, we list the tested tissues with the permeabilization condition used:

=== "Human"

    | Tissue type      | Permeabilization condition (time, pepsin concentration)|
    | ----------- | ----------- | 
    | Metastatic lymph node     | 45 min, 1.4 U/uL| 
    | Healthy  lymph node | 45 min, 1.4 U/uL|
    | Head and neck squamous cell carcinoma| 45 min, 1.4 U/uL | 
    | iPSC-derived Brain Organoids | 15 - 30 min, 0.7 U/uL|

=== "Mouse"

    | Tissue type      | Permeabilization conditions|
    | ----------- | ----------- | 
    | Mouse head (embryo E13)| 30 min, 0.7 U/uL|
    | Mouse brain| 30 min, 0.7 U/uL|

## Capture area
[__Can you store capture area pieces?__](#store-capture-areas){ #store-capture-areas }

Yes. We have successfully generated libraries from capture areas stored >12 months. We recommend storage at -20°C or -80°C with silica beads. 

[__Can tissue sections be placed on a capture area and stored before library preparation?__](#place-store-capture-areas){ #place-store-capture-areas }

Yes. We have stored capture areas with tissue sections for 1 week at -80°C before proceeding with library
preparation and have not observed an effect mRNA capture (qPCR). 

Sections on capture areas were stored before methanol fixation. Once the tissue was placed on the capture area,
care was taken to keep it frozen until library preparation.

Longer storage may be possible, but has not been systematically tested. 

## General protocol
[__What instruments/equipment is required to perform the open-ST protocol?__](#instruments-required){ #instruments-required }

Open-ST was developed with the idea to make it accessible to any laboratory. It requires [standard lab equipment](experimental/getting_started.md): 

-	***Capture area generation***: Illumina® NovaSeq6000, 3D-printer (*capture areas can be made without a 3D-printed cutting guide, if this is not available. A cutting guide is recommended for ease-of-use, preventing scratches and irregular capture areas*)  
- 	***Library preparation***: Cryostat, hybridization oven, thermocycler, Bluepippin or PippinHT (alternatively, agarose gel and DNA extraction can be done manually)
-	***Imaging***: Brightfield microscope / Fluorescence microscope 
-	***Quality control***: Qubit, qPCR, automated gel electrophoresis (Bioanalyzer, TapeStation, Fragment analyzer)
- 	***Sequencing***: various sequencers can be used, as long as a minimum of 100 cylces can be sequenced (ex., Illumina® NextSeq500/550, NextSeq2000, NovaSeq6000, NovaSeqX(plus))

## Imaging
[__Are immunohistochemical (IHC) or immunofluorescence (IF) stainings compatible with open-ST?__](#imaging-ihc-if){ #imaging-ihc-if }

IHC/IF staining may reduce the quality of the resulting open-ST library, since staining occurs before RNA capture and may lead to RNA degradation.  

We have successfully applied hematoxilin and eosin (H&E) staining, as well as fluorescent cytoskeletal (Phalloidin) and nuclei (DAPI) staining, as part of the open-ST workflow.


## Library preparation and sequencing
[__How does a good library profile look like before and after size selection?__](#what-good-profile){ #what-good-profile }
		
A good library profile is smooth without any short fragment peaks or evident peaks inside the library range.
If small length peaks remain after library size selection, they should be removed with an additional size selection (agarose gel or beads). 

Evident peaks inside the library range could be due to over-amplification of a low complexity library.
If possible, we recommend checking previous protocol steps, including RNA quality control, permeabilization condition, amplification cycling number.

[__What is the recommended sequencing depth?__](#what-sequencing-depth){ #what-sequencing-depth }

Sequencing depth requirements vary with tissue section size, coverage of capture area and experimental goals.
For reference, for a 3x4 mm section we obtain a median of around 900 UMIs per cell when investing 500 million sequencing reads.  
Shallow sequencing can always be performed first to assess general library quality (%spatial mapping, % uniquely mapping to genome, % rRNA, % mt-encoded). 


## Image processing
[__Alignment of H&E and spatial (visibility of alignment marks, how to ensure, how many necessary?)__](#fiducial-alignment){ #fiducial-alignment }

TODO: @Daniel

## Sequencing processing
[__Perhaps smth on index hopping__](#index-hopping){ #index-hopping }

TODO: @Daniel
https://www.10xgenomics.com/blog/answering-your-questions-about-the-visium-spatial-gene-expression-solution 