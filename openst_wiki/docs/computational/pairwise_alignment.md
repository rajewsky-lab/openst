# Pairwise alignment of imaging and spatial transcriptomics
In order to assign transcripts to the nuclei segmented from the staining images, the pairwise alignment between the imaging and spatial transcriptomics modality must be performed. For this, we provide software that allows to:

1. Create pseudoimages from the spatial transcriptomics data.
2. Align H&E images to pseudoimages of ST data in two steps (coarse, fine).
3. Manually validate and refine the alignment for more precision, with an interactive user interface.

## Pairwise alignment
```bash
openst pairwise_aligner --args --metadata=<where_to_write_metadata>
```
## Quality Control
```bash
openst report --metadata=<where_to_write_metadata> --output=<path_to_html_file>
```
## (Optional) Refinement of alignment
Use the provided notebook that spawns a napari environment for visual assessment and manual refinement of the alignment, if necessary.
```bash
openst alignment_refiner --args --metadata=<where_to_write_metadata>
```

We can create the report again, to finally validate the refinement as we did for the pairwise alignment
```bash
openst report --metadata=<where_to_write_metadata> --output=<path_to_html_file>