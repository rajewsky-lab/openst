

Troubleshooting

# FAQs 

#### 1. Possibilities for multiplexing (diff tissues, perm times, …)? 
	
Multiplexing your experiments is possible with the use of an Array Hybridization chamber cassette (https://arrayjet.com/consumables/gel-company-arrayslide-16-hybridization-chamber-cassette). You can use up to 16 FC pieces per cassette, that means that you can test in one experiment different condition (i.e. permeabilization time/ concentration) or different sample types.
This will allow you to save experimental time while giving higher reliability in sample/condition comparison. Moreover, open-ST allow you to multiplexing sample with also in NGS run using different indexing during library preparation (_see protocol at step X_@Marie)

#### 2. What tissues have already been tested with open-ST? 

Several tissues have been tested using open-ST, including in vitro 3D cultures. Here we list the tested tissues with the permeabilization condition used:

> Human 

| Tissue type      | Permeabilization conditions|
| ----------- | ----------- | 
| Metastatic lymph node     | 45 min, 1.4 U/uL 
| Metastatic  lymph node | 45 min, 1.4 U/uL
| Lung cancer  | X  min, 0.7 U/uL @Marie
| iPs Brain Organoids | 15 - 30 min, 0.7 U/uL

> Mouse model

| Tissue type      | Permeabilization conditions|
| ----------- | ----------- | 
| Mouse head (embryo E13)| 30 min, 0.7 U/uL
| Mouse brain| 30 min, 0.7 U/uL
| Brest cancer model Stage 1-4 | 30 min - 45 min, 0.7 U/uL @Marie

#### 3. What tissue section thickness is recommended? @Marie

For all experiment we have done we used and thus we recommend to use 10 um tissue thickness. It is possible to try also less than 10 mm or to go up to 12 mm for particular tissue (such as brain). 
It is of crucial importance that the tissue structure is kept intact, thus avoid any tissue rip or fold.

#### 4. Can I do IHC staining on my tissue on the capture area? 

Perform a proper IHC staining or Immunofluorescence because it required several steps (including an overnight incubation) that could degraded RNA. 
So it will be challenging to implement those protocol on open-ST workflow. 
However, we have been testing staining that required only few steps such as Hematoxylin and Eosin to retrieve tissue morphology or a fluorescent based staining to define cell membrane (Phalloidin) and nuclei (DAPI/Hoechst) staining.

#### 5. Do you recommend replicates? Biological/ technical 

Yes, is important to have replicate expecially biological due to the intervariability between samples. However, a tecnical replicate could be avoided due to the high reproducibility of our protocol.

#### 6. What is the recommended freezing and embedding process? 

An optimal embedding process is required to avoid freezing artifact, that are frequent especially in less dense tissues. The best results are achieved when tissue is frozen in a fast and uniform way.
Thus, the recommended procedure consist in simultaneous freezing and OCT-embedded of the tissue using an isopentane bath in dry ice or liquid nitrogen. 
However, is also possible to embed already snap-frozen tissue in OCT or use of isopentane using power of dry ice. In this cases you should be aware of:

- Do not immerge the tissue directly on liquid nitrogen or dry ice (it will burn the edges!)
- Once tissue is frozen avoid to melt the tissue wen embedding in OCT
- Keep the embedding mold bases fully covered by dry ice to allow an homogeneous freezing

#### 7. What instrument do I need in lab to perform open-ST protocol?

Open-ST was built with the idea to make it accessible to any laboratory, so there are no specific with the use of routine lab instrumentation.

-	***Library preparation***: Hybridization oven system, Thermocycler, Bluepippin
-	***Imaging***: Brightfield microscope / Fluorescence microscope 
-	***Quality control***: qPCR and TapeStationn

#### 8. Can a monochrome camera be used? 

A monochrome camera can be used only on fluorescent microscopy if it is used the immunofluorescent protocol to retrieve tissue morphology (i.e. Phalloidin and Dapi). However, for H&E protocol is required a color camera (if feasible with good resolution, scanning and with balancing functionality)
Imaging https://kb.10xgenomics.com/hc/en-us/articles/360035999152-What-are-the-imaging-system-requirements-for-running-Visium-for-fresh-frozen- 
		
#### 9. What pepsin incubation times should I test? 

Is important to set the pepsin incubation time and concentration for each tissue types. 
We recommend to test at least a minimum range including 15 min, 30 min, and 60 min with two different concentration (0.7 and 1.4 U/uL ). 
It is preferable to chose the minimum incubation times that gives the maximum RNA recover (_see the protocol “qPCR/permeabilization condition”_ @Marie)

#### 10. What do I do if tissue remains after the tissue removal step? 

Normally is all the steps should be removed efficiently (@Marie I don’t remember in which condition the tissue was difficult to remove but was I remember that was due to some specific steps…permeabilization? do you remember? I look through my notes) .
However if you still see trace of tissue after removal step you can add additional washes to remove it completely.

#### 11.How does a good library profile look like before and after size selection? 
		
A good library profile is smooth without any short peaks or evident peaks inside the library range. If small length peaks remain after library size selection, it is important to remove with an additional size selection (BuePippin or beads). Evident peak inside the library range could be due to over amplification artifact, if possible we recommend to check all the previous steps, permeabilization condition , amplification cycling number.

**Figure of good library \Figure of over amplificated library and before Bluepippin**

#### 12.What are the best practices to avoid cross-contamination? 

In order to avoid any RNA contamination is important to wipe down cryostat, brushes with stage with 80%-100% ethanol and change blade for different tissues.

#### 13.Can tissue sections be placed on a capture area and stored before library preparation? 

Yes, you can store tissue placed on capture area for up to @Marie months. On dry or methanol? Before or after methanol fixation?

#### 14.Is it possible to store capture area pieces?

Yes, you can store it dry at -80 with silica beads for up to 6 months (?@Marie)

#### 15.What is the recommended sequencing depth? 

-	https://kb.10xgenomics.com/hc/en-us/articles/360036045332-What-are-the-minimum-recommended-sequencing-specifications-for-Visium-for-fresh-frozen-libraries- 
-	QC through shallow seq? What QC metrics are valid at low seq depth? https://kb.10xgenomics.com/hc/en-us/articles/12769731844621-Can-I-perform-shallow-sequencing-to-assess-the-quality-of-Visium-Spatial-Gene-Expression-for-Fresh-Frozen-libraries- 
			
#### 16.Can open-ST libraries be sequenced together with other sequencing libraries? 

@Marie

#### 17.Alignment of H&E and spatial (visibility of alignment marks, how to ensure, how many necessary?) 

@Daniel



Perhaps smth on index hopping

https://www.10xgenomics.com/blog/answering-your-questions-about-the-visium-spatial-gene-expression-solution 


