import numpy as np
import os

from PIL import Image
from tqdm import tqdm

# from openst.preprocessing.CUT.models import create_model

def _images_to_tiles(img):
    _img_shape = img.shape
    tiles = []
    for x in range(0, _img_shape[0]-args.tile_size_px, args.tile_size_px):
        for y in range(0, _img_shape[1]-args.tile_size_px, args.tile_size_px):
            tiles.append([x, y])

    imgs = []
    for i, coord in tqdm(enumerate(tiles)):
        imgs.append(img[coord[0]:(coord[0]+args.tile_size_px), coord[1]:(coord[1]+args.tile_size_px)])

    return tiles, imgs 


def _tiles_to_images(tiles, imgs, dest_shape):
    img_restitch = np.zeros(dest_shape)
    for i, coord in tqdm(enumerate(tiles)):
        img_restitch[coord[0]:(coord[0]+args.tile_size_px), coord[1]:(coord[1]+args.tile_size_px)] = imgs[i] 


def _image_preprocess(imgs_tiles, model):
    for i, data in enumerate(dataset):
        if i == 0:
            model.data_dependent_initialize(data)
            model.setup(opt)
            model.parallelize()
            if args.eval:
                model.eval()
        data = {'A': None, 'B': None, 'A_paths': None, 'B_paths': None}
        model.set_input(data)  # unpack data from data loader
        model.test()           # run inference
        visuals = model.get_current_visuals()  # get image results


def _run_image_preprocess():
    # loading arguments
    args = parser.parse_args()

    Image.MAX_IMAGE_PIXELS = 933120000

    # Load image and get shape
    img = np.array(Image.open(args.input_img))
    _img_shape = img.shape

    tiles, imgs_tiles = _images_to_tiles(img)

    opt = TestOptions().parse()
    opt.num_threads = 1
    opt.batch_size = 1
    opt.serial_batches = True
    opt.no_flip = True
    opt.load_size = args.tile_size_px

    model = create_model(opt)
    imgs_tiles_processed = _image_preprocess(imgs_tiles)

    img_restitch = _tiles_to_images(tiles, imgs_tiles_processed, _img_shape)
    Image.fromarray(img_restitch.astype(np.uint8)).save(args.output_img)

if __name__ == "__main__":
    from openst.cli import get_image_preprocess_parser
    args = get_image_preprocess_parser().parse_args()
    _run_image_preprocess(args)