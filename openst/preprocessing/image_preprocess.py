import numpy as np
import h5py
import logging
import pathlib
import os
import torch

from PIL import Image
from tqdm import tqdm

from torch.utils.data import DataLoader
import torchvision.transforms as transforms

from openst.preprocessing.CUT.models import create_model
from openst.preprocessing.CUT.options.test_options import TestOptions
from openst.utils.file import check_directory_exists, check_file_exists, download_url_to_file

OPENST_MODEL_NAMES = [
    "HE_CUT_rajewsky"
]

MODEL_DIR = pathlib.Path.home().joinpath(".CUT", "models")
_MODEL_URL = "http://bimsbstatic.mdc-berlin.de/rajewsky/openst-public-data/CUT_models"

def create_dataset(opt):
    dataset = OpenSTDataset(opt)
    dataloader = DataLoader(
            dataset,
            batch_size=1,
            shuffle=False,
            num_workers=1,
            drop_last=False,
        )
    return dataset

def cache_model_path(basename):
    import tarfile

    if basename not in OPENST_MODEL_NAMES:
        logging.error(f"The model {basename} is not available. It has to be one of {OPENST_MODEL_NAMES}")
        exit(1)

    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    url = f"{_MODEL_URL}/{basename}.tar.xz"
    cached_file = os.fspath(MODEL_DIR.joinpath(basename))
    if not os.path.exists(cached_file):
        logging.info('Downloading: "{}" to {}'.format(url, cached_file))
        download_url_to_file(url, f"{cached_file}.tar.xz", progress=True)

        logging.info('Decompressing "{}" to {}'.format(f"{cached_file}.tar.xz", MODEL_DIR))
        with tarfile.open(f"{cached_file}.tar.xz", "r:xz") as tar:
            tar.extractall(path=MODEL_DIR)

    return cached_file

def get_transform():
    transform_list = []
    transform_list += [transforms.ToTensor()]
    transform_list += [transforms.Normalize((0.5, 0.5, 0.5), (0.5, 0.5, 0.5))]

    return transforms.Compose(transform_list)

class OpenSTDataset():
    """Wrapper class of Dataset class that performs multi-threaded data loading"""

    def __init__(self, image_tiles):
        self.image_tiles = image_tiles

    def __len__(self):
        """Return the number of data in the dataset"""
        return len(self.image_tiles)
    
    def __getitem__(self, index):
        """Return a data point and its metadata information.

        Parameters:
            index (int)      -- a random integer for data indexing

        Returns a dictionary that contains A, B, A_paths and B_paths
            A (tensor)       -- an image in the input domain
            B (tensor)       -- its corresponding image in the target domain
            A_paths (str)    -- image paths
            B_paths (str)    -- image paths
        """
        transform = get_transform()
        A = transform(self.image_tiles[index])
        B = None

        return {'A': A[None], 'B': B, 'size': A[None].shape}

def _images_to_tiles(img, tile_size_px=512):
    _img_shape = img.shape
    tiles = []
    for x in range(0, _img_shape[0]-tile_size_px, tile_size_px):
        for y in range(0, _img_shape[1]-tile_size_px, tile_size_px):
            tiles.append([x, y])

    imgs = []
    for i, coord in tqdm(enumerate(tiles)):
        imgs.append(img[coord[0]:(coord[0]+tile_size_px), coord[1]:(coord[1]+tile_size_px)])

    return tiles, imgs 


def _tiles_to_images(tiles, imgs, dest_shape, tile_size_px=512):
    img_restitch = np.zeros(dest_shape)
    for i, coord in tqdm(enumerate(tiles)):
        img_restitch[coord[0]:(coord[0]+tile_size_px), coord[1]:(coord[1]+tile_size_px)] = imgs[i] 


def _image_preprocess(model, dataset):
    model.data_dependent_initialize(dataset[0])
    model.eval()
    
    output = []
    for data in tqdm(dataset):
        model.set_input(data)
        model.test()
        output += [model.get_current_visuals()]
    
    return output

def _load_image_adata(args):
    if args.h5_in != '':
        check_file_exists(args.h5_in)
        adata = h5py.File(args.h5_in, 'r+')
        im = adata[args.image_in][:]
    else:
        check_file_exists(args.image_in)
        if not check_directory_exists(args.mask_out, exception=False):
            raise FileNotFoundError("Parent directory for --mask-out does not exist")
        
        Image.MAX_IMAGE_PIXELS = 933120000
        im = np.array(Image.open(args.image_in))        
        adata = None
        
    return adata, im

def _save_image_adata(adata, restored_img, args):
    if args.h5_in != "":
        if args.image_out in adata:
            logging.warn(f"The object {args.image_out} will be removed from the h5py file")
            del adata[args.image_out]

        logging.info(f'Saving mask to adata in {args.image_out}')
        adata[args.image_out] = restored_img
    else:
        logging.info(f'Saving mask to separate file in {args.image_out}')
        Image.fromarray(restored_img).save(args.image_out)

def _run_image_preprocess(args):
    if args.h5_in == "" and args.image_in == "":
        logging.error("You need to provide at least one of `--h5-in` or `--image-in`")
        exit(1)

    adata, img =_load_image_adata(args)
    _img_shape = img.shape

    logging.info(f"Converting image of shape {_img_shape} into tiles")
    tiles, imgs_tiles = _images_to_tiles(img, args.tile_size_px)

    opt = TestOptions("").parse()
    opt.num_threads = 1
    opt.batch_size = 1
    opt.serial_batches = True
    opt.no_flip = True
    opt.load_size = args.tile_size_px
    opt.crop_size = args.tile_size_px
    opt.pretrained_name = None

    model = create_model(opt)
    model.save_dir = cache_model_path(args.model)
    model.setup(opt)
    model.parallelize()

    logging.info(f"Creating dataset from image ({len(imgs_tiles)} tiles)")
    dataset = create_dataset(imgs_tiles)
    
    logging.info(f"Running restoration on {len(imgs_tiles)} tiles")
    imgs_tiles_processed = _image_preprocess(model, dataset)

    logging.info("Merging back into single image")
    restored_img = _tiles_to_images(tiles, imgs_tiles_processed, _img_shape, args.tile_size_px)
    _save_image_adata(adata, restored_img, args)

if __name__ == "__main__":
    from openst.cli import get_image_preprocess_parser
    args = get_image_preprocess_parser().parse_args()
    _run_image_preprocess(args)