import h5py
import logging
import dask.array as da
import dask_image.imread

def _run_merge_modalities(args):
    logging.info(f"Updating {args.h5_in} in place")
    image = dask_image.imread.imread(args.image_in)[0]

    with h5py.File(args.h5_in, 'r+') as adata:
        if args.image_key in adata:
            logging.warn(f"The key '{args.image_key}' will be overwritten")
            del adata[args.image_key]

        dset = adata.create_dataset(args.image_key, shape=image.shape,
                                        dtype=image.dtype)  
        logging.info(f"Saving mask to Open-ST h5 object in key '{args.image_key}'")
        da.store(image, dset)

if __name__ == "__main__":
    from openst.cli import get_merge_modalities_parser
    args = get_merge_modalities_parser().parse_args()
    _run_merge_modalities()


    