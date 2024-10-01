import numpy as np
import h5py
import logging
import pathlib
import os

from PIL import Image
from tqdm import tqdm
from threadpoolctl import threadpool_limits

try:
    import torchvision.transforms as transforms
except ImportError:
    raise ImportError(
        "Please install napari: `pip install torchvision` "+
        "or find more information at https://pytorch.org/get-started/locally/"
    )

from openst.preprocessing.CUT.models import create_model
from openst.preprocessing.CUT.options.test_options import TestOptions
from openst.preprocessing.CUT.util import tensor2im
from openst.utils.file import check_directory_exists, check_file_exists, download_url_to_file

OPENST_MODEL_NAMES = [
    "HE_CUT_rajewsky"
]

MODEL_DIR = pathlib.Path.home().joinpath(".openst", "CUT", "models")
_MODEL_URL = "http://bimsbstatic.mdc-berlin.de/rajewsky/openst-public-data/CUT_models"

def create_dataset(args, opt):
    dataset = OpenSTDataset(args, opt)
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
            tar.extractall(path=MODEL_DIR, filter='data')

    return cached_file

def get_transform(args):
    transform_list = []

    osize = [args.tile_size_px]*2
    transform_list.append(transforms.Resize(osize, Image.BICUBIC))
    transform_list.append(transforms.RandomCrop(args.tile_size_px))
    transform_list.append(transforms.Lambda(lambda img: __make_power_2(img, base=4, method=Image.BICUBIC)))
    transform_list += [transforms.ToTensor()]
    transform_list += [transforms.Normalize((0.5, 0.5, 0.5), (0.5, 0.5, 0.5))]

    return transforms.Compose(transform_list)

def __make_power_2(img, base, method=Image.BICUBIC):
    ow, oh = img.size
    h = int(round(oh / base) * base)
    w = int(round(ow / base) * base)
    if h == oh and w == ow:
        return img

    return img.resize((w, h), method)

class OpenSTDataset():
    """Wrapper class of Dataset class that performs multi-threaded data loading"""

    def __init__(self, args, image_tiles):
        self.image_tiles = image_tiles
        self.args = args

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
        transform = get_transform(self.args)
        A = transform(self.image_tiles[index])
        B = None

        return {'A': A[None], 'B': B, 'size': A[None].shape}

def _image_to_tiles(img, tile_size_px=512):
    _img_shape = img.shape
    tiles = []
    for x in range(0, _img_shape[0]-tile_size_px, tile_size_px):
        for y in range(0, _img_shape[1]-tile_size_px, tile_size_px):
            tiles.append([x, y])

    imgs = []
    for coord in tqdm(tiles):
        imgs.append(Image.fromarray(img[coord[0]:(coord[0]+tile_size_px), coord[1]:(coord[1]+tile_size_px)]).convert('RGB'))

    return tiles, imgs 


def _tiles_to_image(tiles, imgs, dest_shape, tile_size_px=512):
    img_restitch = np.zeros(dest_shape, dtype=np.uint8)
    for i, coord in tqdm(enumerate(tiles)):
        img_restitch[coord[0]:(coord[0]+tile_size_px), coord[1]:(coord[1]+tile_size_px)] = imgs[i] 

    return img_restitch


def _image_preprocess(model, dataset):
    model.data_dependent_initialize(dataset[0])
    model.eval()
    
    output = []
    for data in tqdm(dataset):
        model.set_input(data)
        model.test()
        # we need to transform (C, X, Y) to original (X, Y, C)
        _tensor_output = model.get_current_visuals()['fake_B'].cpu()
        _image_output = tensor2im(_tensor_output)
        output += [_image_output]
    
    return output

def _load_image_adata(args):
    if args.h5_in != '':
        check_file_exists(args.h5_in)
        adata = h5py.File(args.h5_in, 'r+')
        im = adata[args.image_in][:]
    else:
        check_file_exists(args.image_in)
        if not check_directory_exists(args.image_out, exception=False):
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

def run_image_preprocess(args):
    if args.h5_in == "" and args.image_in == "":
        logging.error("You need to provide at least one of `--h5-in` or `--image-in`")
        exit(1)

    adata, img = _load_image_adata(args)
    _img_shape = img.shape

    gpu_id = "0" if args.device == "cuda" else "-1"
    opt = TestOptions(f"--gpu_ids {gpu_id}").parse()
    opt.num_threads = 1
    opt.batch_size = 1
    opt.serial_batches = True
    opt.no_flip = True
    opt.load_size = args.tile_size_px
    opt.crop_size = args.tile_size_px
    opt.pretrained_name = None

    logging.info(f"Loading model {args.model}")
    model = create_model(opt)
    model.save_dir = cache_model_path(args.model)
    model.setup(opt)
    model.parallelize()

    logging.info(f"Converting image of shape {_img_shape} into tiles")
    tiles, imgs_tiles = _image_to_tiles(img, args.tile_size_px)

    logging.info(f"Creating dataset from image ({len(imgs_tiles)} tiles)")
    dataset = create_dataset(args, imgs_tiles)
    
    logging.info(f"Running restoration on {len(imgs_tiles)} tiles")
    imgs_tiles_processed = _image_preprocess(model, dataset)

    logging.info(f"Merging {len(tiles)} tiles back into single image of shape {_img_shape}")
    restored_img = _tiles_to_image(tiles, imgs_tiles_processed, _img_shape, args.tile_size_px)
    _save_image_adata(adata, restored_img, args)

def _run_image_preprocess(args):
    with threadpool_limits(limits=args.num_workers):
        run_image_preprocess(args)

if __name__ == "__main__":
    from openst.cli import get_image_preprocess_parser
    args = get_image_preprocess_parser().parse_args()
    _run_image_preprocess(args)