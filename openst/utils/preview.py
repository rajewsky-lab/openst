def _run_preview(args):
    import h5py
    try:
        import napari
    except ImportError:
        raise ImportError(
            "Please install napari: `pip install napari`."
        )

    from openst.utils.file import check_file_exists
    
    check_file_exists(args.adata)
    adata = h5py.File(args.adata, 'r+')

    viewer = napari.Viewer()

    if args.image_keys is not None:
        for i, _image_key in enumerate(args.image_keys):
            if len(args.image_resampling) != len(args.image_keys):
                image_resampling = args.image_resampling[0]
            else:
                image_resampling = args.image_resampling[i]

            viewer.add_image(data=adata[_image_key][::image_resampling, ::image_resampling])

    if args.spatial_coord_keys is not None:
        for i, _spatial_coord_key in enumerate(args.spatial_coord_keys):
            if len(args.spatial_coord_resampling) != len(args.spatial_coord_keys):
                spatial_coord_resampling = args.spatial_coord_resampling[0]
            else:
                spatial_coord_resampling = args.spatial_coord_resampling[i]

            viewer.add_points(data=adata[_spatial_coord_key][::spatial_coord_resampling])
    
    napari.run()

if __name__ == "__main__":
    from openst.cli import get_preview_parser
    args = get_preview_parser().parse_args()
    _run_preview()