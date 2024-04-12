# Align image to transcriptome

In the previous step, the transcriptomic reads were processed and mapped in tissue space with `spacemake`.
Now, we perform a pairwise alignment between the imaging and spatial transcriptomics modality, such that
we can later aggregate transcripts into individual cells delimited by the segmentation mask.

We will illustrate how to do this in a semiautomatic manner: that is, running the coarse alignment in
an automatic fashion, and the fine alignment (to fiducial marks) via GUI, in a manual manner. Although we
provide models for fiducial feature detection, the accuracy might be affected by the type of microscope,
imaging strategy, tissue type and width... Thus, manual fine alignment is a good option. Thanks to the
GUI specifically designed for this task, the necessary amount of time invested is reasonable
(~5 minutes per sample consisting of 12 tiles).

## Coarse alignment
We assume that you have the single `h5ad` file that contains all the tiles for the experiment (you got it automatically
from `spacemake`, or you followed the [stitching instructions](preprocessing_sequencing.md#expected-output)).

Now, make sure of two things:
1. You have an environment where the `openst` package is installed (we recommend doing this
in an environment different to the one used for `spacemake`).
1. You are in the `/home/user/openst_e13_demo` folder (or similar, depending on your system, username...)

```sh
openst pairwise_aligner \
    --image-in alignment/e13_mouse_head.tif \
    --h5-in projects/openst_demo/processed_data/openst_demo_e13_mouse_head/illumina/complete_data/dge/dge.all.polyA_adapter_trimmed.mm_included.spatial_beads_puck_collection.h5ad \
    --h5-out alignment/openst_demo_e13_mouse_head_spatial_beads_puck_collection_aligned.h5ad \
    --save-image-in-h5 \
    --only-coarse \
    --device cuda
```

This will create the file `alignment/openst_demo_e13_mouse_head_spatial_beads_puck_collection_aligned.h5ad`, that hopefully contains a
proper coarse alignment of the two modalities (at low resolution).

## Fine alignment (manual)

You will now visually asses it and refine the previous coarse alignment using the provided GUI tool. 
From the terminal, run:

```sh
openst manual_pairwise_aligner_gui
```

Then, follow these instructions on the GUI, by pressing the *buttons*:

1. Click on *Load h5ad* and browse to the location where the `openst_demo_e13_mouse_head_spatial_beads_puck_collection_aligned.h5ad` file
is located. Select it and click *Open*. 
2. Go to *Data properties*, then click under *Image data*, and select from the tree `uns>spatial_pairwise_aligned>staining_image_transformed`. Click *Ok*.
3. Click under *Spatial coordinates*, and select from the tree `obsm/spatial_pairwise_aligned_coarse`. Click *Ok*.
4. Go to rendering settings, and keep in mind the *Threshold counts* parameter.
5. Now select the *all_tiles_coarse* option from the *Layer selector*, and click on *Render*. You will see the staining image on the upper left,
   the transcriptomic image on the top right, and the merge on the bottom left. These two modalities should *roughly* match. If not, you would need to
   run the coarse alignment in manual mode, too.
6. Once you've checked that the coarse alignment is fine, start selecting the layer 0, and click *Render*.
7. Now, you need to select the pairs of corresponding fiducial markers on both modalities (top left and top right) by double clicking on the left image, first, and on the right image, second. You can drag the points to new locations, zoom into the images with the mouse wheel, pan the image by holding the mouse right cursor and moving, and remove the points by pressing backspace on your keyboard.
8. You need to select at least 2 pairs of points per layer. Don't forget to click on *Render* after you select a new layer (e.g., 1, 2, 3...) so it displays. Also, you can preview how the alignment will look like after transformation by pressing *Preview alignment*. This will update the merge view (bottom left). You can change the opacity of the merge images under *Rendering settings*, by adjusting the Image A/B opacity sliders. Also, if you don't see a lot of signal in the transcriptomic pseudoimage, try decreasing the *Threshold counts* and *Pseudoimage size* values. Per tile (layer), adjust the keypoints until you get a good alignment as many times as necessary.
9. Once you've finished, go under *Keypoint properties* and click on *Save keypoints*. Go to the `alignment` folder under `openst_e13_demo`, and save the file as `openst_e13_demo_fine_keypoints`.

Great! You have just created a file with the corresponding points between the respective coordinate systems of
the spatial trancriptome and image modalities. 

## Apply keypoint transformation to coarse alignment
Now, you need to run a program that takes the `h5ad` file with the
coarse alignment, and the keypoints file, to perform the fine alignment:

```sh
openst manual_pairwise_aligner \
    --keypoints-in alignment/openst_e13_demo_fine_keypoints.json \
    --h5-in alignment/openst_demo_e13_mouse_head_spatial_beads_puck_collection_aligned.h5ad \
    --fine
```

After this, no file will be created nor removed; the coordinates of the fine alignment will be added to the existing
`openst_demo_e13_mouse_head_spatial_beads_puck_collection_aligned.h5ad` file.

That's it! Now you're ready to go to the next step.