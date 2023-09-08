import matplotlib.pyplot as plt
import os
from multiprocessing import Pool
from tqdm import tqdm
from PIL import Image

import numpy as np
from cellpose import models
from cellpose.io import imread

Image.MAX_IMAGE_PIXELS = 933120000

INPUT_IMG='/data/rajewsky/home/dleonpe/projects/openst_paper/data/1_images_mouse/fc_sts_75_1b/stitched/Image_Stitched_Composite.tif'
OUTPUT_MASK = '/data/rajewsky/home/dleonpe/projects/openst_paper/data/1_images_mouse/fc_sts_75_1b/segmented/Image_Stitched_Composite_mask.tif'
MODEL_PATH = '/data/rajewsky/home/dleonpe/projects/openst_paper/data/0_cellpose_models/HE_cellpose_rajewsky'
FLOW_THRESHOLD=2
PROB_THRESHOLD=-0.5
SQUARE_SIZE = 1024
OFFSET = 1024
MODEL_DIAMETER = 20

# Load cellpose model
model = models.CellposeModel(gpu=True, pretrained_model=MODEL_PATH)

# Load WSI image
img_wsi = np.array(Image.open(INPUT_IMG))
imgs = []

# Create tiles
for x in range(0, img_wsi.shape[0]-SQUARE_SIZE, OFFSET):
    for y in range(0, img_wsi.shape[1]-SQUARE_SIZE, OFFSET):
        imgs.append(img_wsi[x:(x+SQUARE_SIZE), y:(y+SQUARE_SIZE)])

# Run prediction per tile
masks = []
for im in tqdm(imgs):
    channels = [[0,0]]
    mask, _, _ = model.eval([im], diameter=MODEL_DIAMETER, channels=channels, flow_threshold=2, cellprob_threshold=-0.5)
    masks.append(mask[0])

# Stitch mask back into WSI image
img_restitch = np.zeros(img_wsi.shape[0:2])

_image_index = 0
for x in range(0, img_wsi.shape[0]-SQUARE_SIZE, OFFSET):
    for y in range(0, img_wsi.shape[1]-SQUARE_SIZE, OFFSET):
        img_restitch[x[0]:(x[0]+SQUARE_SIZE), y[1]:(y[1]+SQUARE_SIZE)] = masks[_image_index]
        _image_index += 1

Image.fromarray(img_restitch.astype(np.uint32)).save(OUTPUT_MASK)