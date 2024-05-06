import h5py
import logging

from openst.utils.pseudoimage import create_unpaired_pseudoimage
from openst.utils.file import h5_to_dict

def _show_viewer(adata, args):
    try:
        import napari
    except ImportError:
        raise ImportError(
            "Please install napari: `pip install napari`."
        )

    viewer = napari.Viewer()

    if args.image_keys is not None:
        for i, _image_key in enumerate(args.image_keys):
            if len(args.image_resampling) != len(args.image_keys):
                image_resampling = args.image_resampling[0]
            else:
                image_resampling = args.image_resampling[i]

            viewer.add_image(data=adata[_image_key][::image_resampling, ::image_resampling])

    if args.pseudoimage_keys is not None:
        for i, _pseudoimage_key in enumerate(args.pseudoimage_keys):
            if len(args.pseudoimage_units_to_um) != len(args.pseudoimage_keys):
                pseudoimage_units_to_um = args.pseudoimage_units_to_um[0]
            else:
                pseudoimage_units_to_um = args.pseudoimage_units_to_um[i]

            _pseudoimage, _ = create_unpaired_pseudoimage(adata, _pseudoimage_key, pseudoimage_units_to_um, write_rescaled=False)

            viewer.add_image(data=_pseudoimage)

    if args.spatial_coord_keys is not None:
        for i, _spatial_coord_key in enumerate(args.spatial_coord_keys):
            if len(args.spatial_coord_resampling) != len(args.spatial_coord_keys):
                spatial_coord_resampling = args.spatial_coord_resampling[0]
            else:
                spatial_coord_resampling = args.spatial_coord_resampling[i]

            viewer.add_points(data=adata[_spatial_coord_key][::spatial_coord_resampling])
    
    napari.run()

def _show_tree(adata, args):
    import json

    adata_tree = h5_to_dict(adata)
    d = json.dumps(adata_tree, indent=4)
    print(d)
    

def _run_preview(args):
    from openst.utils.file import check_file_exists
    
    check_file_exists(args.h5_in)
    adata = h5py.File(args.h5_in, 'r')

    if args.file_structure:
        logging.info(f"Showing structure from {args.h5_in}")
        _show_tree(adata, args)
    else:
        logging.info(f"Opening napari for {args.h5_in}")
        _show_viewer(adata, args)

if __name__ == "__main__":
    from openst.cli import get_preview_parser
    args = get_preview_parser().parse_args()
    _run_preview()