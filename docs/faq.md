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

The optimal cutting temperature depends on the tissue being sectioned. As a starting point, page 91 in this [reference](../387928-EP_05 - NX70 Op Guide.pdf) is helpful.

## Sample handling

[__Are technical replicates required?__](#technical-replicates){ #technical-replicates }

As we have observed a high reproducibility across technical replicates (consecutive sections), technical replicates are not essential.

[__How many sections (conditions/sample types) can be processed at once?__](#samples-process-at-once){ #samples-process-at-once }
	
Multiple samples can be processed at once with the use of a multiwell chamber hybridization cassette
(such as this one by [ArrayJet](https://arrayjet.com/consumables/gel-company-arrayslide-16-hybridization-chamber-cassette)). 

This cassette has sixteen 7 x 7 mm wells, each fitting one capture area; thus, allowing 16 samples or conditions to be processed simultaneously.
We recommend handling a maximum of 15 capture areas per person for protocol steps 3.1 to 3.5 (until reverse transcription overnight incubation).

Moreover, open-ST libraries are indexed on the p7-adapter side, allowing multiplexing of samples within one NGS run. 

[__What are the best practices to avoid cross-contamination?__](#best-practices-x-contamination){ #best-practices-x-contamination }

In order to avoid any RNA contamination is important to wipe down the cryostat and any tools used (brushes, tweezers) before sectioning with 80%-100% ethanol. Additionaly, change the blade and wipe down the cryostat and tools in between sectioning different samples. 

[__What pepsin incubation times should I test?__](#pepsin-incubation){ #pepsin-incubation }

It is important to set the permeabilization condition for each tissue type. 
We recommend to test at least a minimum range including 15 min, 30 min, and 60 min with two different concentration (0.7 and 1.4 U/uL). 

It is preferable to chose the minimum incubation time/ enzyme concentration that
gives the maximum RNA capture (see [Permeabilization](experimental/library_preparation.md#permeabilization)).

[__What tissues have already been tested with open-ST?__](#tested-tissues){ #tested-tissues }

Several tissues have been tested using open-ST, including *in-vitro* 3D-cultures.
Here, we list the tested tissues with the permeabilization condition used:
**should we give details for other tissues (collaborations)?**

=== "Human"

    | Tissue type      | Permeabilization condition (time, pepsin concentration)|
    | ----------- | ----------- | 
    | Metastatic lymph node     | 45 min, 1.4 U/uL| 
    | Healthy  lymph node | 45 min, 1.4 U/uL|
    | Head and neck squamous cell carcinoma| | 
    | Non-small cell lung cancer  | X  min, 0.7 U/uL @Marie|
    | iPSC-derived Brain Organoids | 15 - 30 min, 0.7 U/uL|

=== "Mouse"

    | Tissue type      | Permeabilization conditions|
    | ----------- | ----------- | 
    | Mouse head (embryo E13)| 30 min, 0.7 U/uL|
    | Mouse brain| 30 min, 0.7 U/uL|
    | Breast cancer model Stage 1-4 | 30 min - 45 min, 0.7 U/uL @Marie|

## Capture area
[__Can you store capture area pieces?__](#store-capture-areas){ #store-capture-areas }

Yes. We store them dry at -80°C with silica beads for up to 6 months (?@Marie)

[__Can tissue sections be placed on a capture area and stored before library preparation?__](#place-store-capture-areas){ #place-store-capture-areas }

Yes. We have stored capture areas with tissue sections for 1 week at -80°C before proceeding with library
preparation and have not observed an effect mRNA capture (qPCR). 

Sections on capture areas were stored before methanol fixation. Once the tissue was placed on the capture area,
care was taken to keep it frozen until library preparation.

Longer storage is possible, but has not been systematically tested. 

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

IHC/IF staining may reduce the quality of the resulting open-ST library. 

not been  Perform a proper IHC staining or Immunofluorescence because it required several steps (including an overnight incubation) that could degraded RNA. 
So it will be challenging to implement those protocol on open-ST workflow.

However, we have been testing staining that required only few steps such as Hematoxylin and Eosin to
retrieve tissue morphology or a fluorescent based staining to define cell membrane (Phalloidin) and nuclei (DAPI/Hoechst) staining.

[__Can a monochrome camera be used?__](#monochrome-camera){ #monochrome-camera }

A monochrome camera can be used only on fluorescent microscopy if it is used the immunofluorescent
protocol to retrieve tissue morphology (i.e. Phalloidin and Dapi).

However, a color camera is required for H&E protocol - if feasible with good resolution, scanning and balancing functionality.

## Library preparation and sequencing
[__How does a good library profile look like before and after size selection?__](#what-good-profile){ #what-good-profile }
		
A good library profile is smooth without any short fragment peaks or evident peaks inside the library range.
If small length peaks remain after library size selection, it is important to remove with an additional size selection (agarose gel or beads). 

Evident peaks inside the library range could be due to over-amplification of a low complexity library.
If possible, we recommend checking previous protocol steps, including RNA quality control, permeabilization condition, amplification cycling number.

[__What is the recommended sequencing depth?__](#what-sequencing-depth){ #what-sequencing-depth }

-	https://kb.10xgenomics.com/hc/en-us/articles/360036045332-What-are-the-minimum-recommended-sequencing-specifications-for-Visium-for-fresh-frozen-libraries- 
-	QC through shallow seq? What QC metrics are valid at low seq depth? https://kb.10xgenomics.com/hc/en-us/articles/12769731844621-Can-I-perform-shallow-sequencing-to-assess-the-quality-of-Visium-Spatial-Gene-Expression-for-Fresh-Frozen-libraries- 


## Image processing
[__Alignment of H&E and spatial (visibility of alignment marks, how to ensure, how many necessary?)__](#fiducial-alignment){ #fiducial-alignment }

TODO: @Daniel

## Sequencing processing
[__Perhaps smth on index hopping__](#index-hopping){ #index-hopping }

TODO: @Daniel
https://www.10xgenomics.com/blog/answering-your-questions-about-the-visium-spatial-gene-expression-solution 