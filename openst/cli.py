import argparse
import logging

DEFAULT_REGEX_TILE_ID = "(L[1-4][a-b]_tile_[1-2][0-7][0-9][0-9])"

PSEUDOIMAGE_HELP = "Generate pseudoimages of Open-ST RNA data and visualize using napari"
def get_pseudoimage_parser():
    parser = argparse.ArgumentParser(
        description=PSEUDOIMAGE_HELP,
        allow_abbrev=False,
        add_help=False,
    )
    # Input
    parser.add_argument(
        "--h5-in",
        type=str,
        required=True,
        help="Necessary to create the pseudoimage",
    )

    # RNA density-based pseudoimage
    parser.add_argument(
        "--spatial-coord-key",
        type=str,
        default="obsm/spatial",
        help="Path to the spatial coordinates inside the spatial object (e.g., 'obsm/spatial')",
    )
    parser.add_argument(
        "--input-resolution",
        type=float,
        default=1,
        help="""Spatial resolution of the input coordinates (retrieved from --spatial-coord-key).
              If it is in microns, leave as 1. If it is in pixels, specify the pixel to micron conversion factor.""",
    )
    parser.add_argument(
        "--render-scale",
        type=float,
        default=2,
        help="Size of bins for computing the binning (in microns). For Open-ST v1, we recommend a value of 2.",
    )
    parser.add_argument(
        "--render-sigma",
        type=float,
        default=1,
        help="Smoothing factor applied to the RNA pseudoimage (higher values lead to smoother images)",
    )
    parser.add_argument(
        "--output-resolution",
        type=float,
        default=0.6,
        help="Final resolution (micron/pixel) for the segmentation mask.",
    )
    return parser


def setup_pseudoimage_parser(parent_parser):
    parser = parent_parser.add_parser(
        "pseudoimage",
        help=PSEUDOIMAGE_HELP,
        parents=[get_pseudoimage_parser()],
    )
    parser.set_defaults(func=cmd_run_pseudoimage_visualizer)

    return parser


def cmd_run_pseudoimage_visualizer(args):
    from openst.utils.pseudoimage import _run_pseudoimage_visualizer

    _run_pseudoimage_visualizer(args)


PREVIEW_HELP = "Preview locations (as points) and images of Open-ST data"
def get_preview_parser():
    parser = argparse.ArgumentParser(
        description=PREVIEW_HELP,
        allow_abbrev=False,
        add_help=False,
    )
    # Input
    parser.add_argument(
        "--h5-in",
        type=str,
        required=True,
        help="Necessary to create the pseudoimage",
    )

    # Check file structure
    parser.add_argument(
        "--file-structure",
        default=False,
        action="store_true",
        help="If set, will not open a visualization screen but will return the tree structure of the h5 file",
    )

    # RNA density-based pseudoimage
    parser.add_argument(
        "--spatial-coord-keys",
        type=str,
        nargs="+",
        default=None,
        help="""Path to the spatial coordinates inside the spatial object (e.g., 'obsm/spatial').
                Can be one or many (separated by space)""",
    )

    # Staining image
    parser.add_argument(
        "--image-keys",
        type=str,
        nargs="+",
        default=None,
        help="""Path to the image to be visualized.
              Can be one or many (separated by space)""",
    )

    parser.add_argument(
        "--pseudoimage-keys",
        type=str,
        nargs="+",
        default=None,
        help="""Path to the spatial coordinates inside the 
              spatial object to visualize as pseudoimage.
              Can be one or many (separated by space)""",
    )

    # Resampling before previsualizing
    parser.add_argument(
        "--spatial-coord-resampling",
        type=int,
        nargs="+",
        default=[1],
        help="""Will load every n-th point. Can be one (same for all spatial-coords)
                or many (1-to-1 mapping to the spatial-coord list)""",
    )
    parser.add_argument(
        "--image-resampling",
        type=int,
        nargs="+",
        default=[1],
        help="""Will load every n-th pixel. Can be one (same for all images)
                or many (1-to-1 mapping to the image list)""",
    )
    parser.add_argument(
        "--pseudoimage-units-to-um",
        type=float,
        nargs="+",
        default=[1.0],
        help="""Conversion factor from spatial units to micron,
        before rendering the pseudoimage. Can be one (same for all images)
                or many (1-to-1 mapping to the image list)""",
    )

    return parser


def setup_preview_parser(parent_parser):
    parser = parent_parser.add_parser(
        "preview",
        help=PREVIEW_HELP,
        parents=[get_preview_parser()],
    )
    parser.set_defaults(func=cmd_run_preview)

    return parser

def cmd_run_preview(args):
    from openst.utils.preview import _run_preview

    _run_preview(args)


MERGE_MODALITIES_HELP = "merge_modalities locations (as points) and images of Open-ST data"
def get_merge_modalities_parser():
    parser = argparse.ArgumentParser(
        description=MERGE_MODALITIES_HELP,
        allow_abbrev=False,
        add_help=False,
    )
    # Input
    parser.add_argument(
        "--h5-in",
        type=str,
        required=True,
        help="Input Open-ST h5 object",
    )

    parser.add_argument(
        "--image-in",
        type=str,
        required=True,
        help="Image that will be loaded and written into the Open-ST h5 object",
    )

    # Staining image
    parser.add_argument(
        "--image-key",
        type=str,
        default="uns/spatial/staining_image",
        help="Key in the Open-ST h5 object where the image will be saved.",
    )

    return parser


def setup_merge_modalities_parser(parent_parser):
    parser = parent_parser.add_parser(
        "merge_modalities",
        help=MERGE_MODALITIES_HELP,
        parents=[get_merge_modalities_parser()],
    )
    parser.set_defaults(func=cmd_run_merge_modalities)

    return parser


def cmd_run_merge_modalities(args):
    from openst.preprocessing.merge_modalities import _run_merge_modalities

    _run_merge_modalities(args)

SEGMENT_HELP = "Image (or pseudoimage)-based segmentation with cellpose and (optional) radial extension"
def get_segment_parser():
    parser = argparse.ArgumentParser(
        description=SEGMENT_HELP,
        allow_abbrev=False,
        add_help=False,
    )

    # Data
    parser.add_argument(
        "--image-in",
        type=str,
        help="""Key in the Open-ST h5 object (when --h5-in is specified)
              or path to the file where the mask will be loaded from""",
    )
    parser.add_argument(
        "--h5-in",
        type=str,
        default="",
        help="If specified, image is loaded from h5 (from key --image-in). Segmentation mask is saved there (to --mask-out)",
    )
    parser.add_argument(
        "--mask-out",
        type=str,
        required=True,
        help="""Key in the Open-ST h5 object (when --h5-in is specified)
              or path to the file where the mask will be written into""",
    )
    parser.add_argument(
        "--rna-segment",
        action="store_true",
        help="""Performs segmentation based on local RNA density pseudoimages from sequencing data,
              instead of using a staining image. 
              This assumes coordinates in microns (can be transformed with --rna-segment-input-resolution)""",
    )

    # Cellpose
    parser.add_argument(
        "--model",
        type=str,
        default="",
        help="""cellpose model - either a path or a valid string to pretrained model.""",
    )
    parser.add_argument(
        "--flow-threshold",
        type=float,
        default=0.5,
        help="cellpose's 'flow_threshold' parameter",
    )
    parser.add_argument(
        "--cellprob-threshold",
        type=float,
        default=0,
        help="cellpose's 'cellprob_threshold' parameter",
    )
    parser.add_argument(
        "--diameter",
        type=float,
        default=20,
        help="cellpose's 'diameter' parameter",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=512,
        help="When prediction of the mask runs in separate chunks, this is the chunk square size (in pixels)",
    )
    parser.add_argument(
        "--chunked",
        action="store_true",
        help="If set, segmentation is computed at non-overlapping chunks of size '--chunk-size'",
    )
    parser.add_argument(
        "--max-image-pixels",
        type=int,
        default=933120000,
        help="Upper bound for number of pixels in the images (prevents exception when opening very large images)",
    )
    parser.add_argument(
        "--device",
        type=str,
        default="cpu",
        choices=["cpu", "cuda"],
        help="Device used to run the segmentation model. Can be ['cpu', 'cuda']",
    )

    # Postprocessing
    parser.add_argument(
        "--dilate-px",
        type=int,
        help="Pixels the outlines of the segmentation mask will be extended",
        required=False,
        default=10,
    )
    parser.add_argument(
        "--outline-px",
        type=int,
        help="Objects will be represented as px-width outlines (only if >0)",
        required=False,
        default=0,
    )

    # Preprocessing
    parser.add_argument(
        "--mask-tissue",
        action="store_true",
        help="Tissue (imaging modality) is masked from the background before segmentation",
    )
    parser.add_argument(
        "--tissue-masking-gaussian-sigma",
        type=int,
        default=5,
        help="The gaussian blur sigma used during the isolation of the tissue on the staining image",
    )
    parser.add_argument(
        "--keep-black-background",
        action="store_true",
        help="Whether to set the background of the imaging modalities to white after tissue masking",
    )

    # RNA density-based segmentation configuration
    parser.add_argument(
        "--rna-segment-spatial-coord-key",
        type=str,
        default="obsm/spatial",
        help="Path to the spatial coordinates inside the spatial object (e.g., 'obsm/spatial')",
    )
    parser.add_argument(
        "--rna-segment-input-resolution",
        type=float,
        default=1,
        help="""Spatial resolution of the input coordinates (retrieved from --rna-segment-spatial-coord-key).
              If it is in microns, leave as 1. If it is in pixels, specify the pixel to micron conversion factor.""",
    )
    parser.add_argument(
        "--rna-segment-render-scale",
        type=float,
        default=2,
        help="Size of bins for computing the binning (in microns). For Open-ST v1, we recommend a value of 2.",
    )
    parser.add_argument(
        "--rna-segment-render-sigma",
        type=float,
        default=1,
        help="Smoothing factor applied to the RNA pseudoimage (higher values lead to smoother images)",
    )
    parser.add_argument(
        "--rna-segment-output-resolution",
        type=float,
        default=0.6,
        help="Final resolution (micron/pixel) for the segmentation mask.",
    )

    # General
    parser.add_argument(
        "--num-workers",
        type=int,
        help="Number of CPU workers when --chunked is specified",
        required=False,
        default=-1,
    )
    parser.add_argument(
        "--metadata",
        type=str,
        default="",
        help="""Path where the metadata will be stored.
        If not specified, metadata is not saved.
        Warning: a report (via openst report) cannot be generated without metadata!""",
    )
    return parser


def setup_segment_parser(parent_parser):
    parser = parent_parser.add_parser(
        "segment",
        help=SEGMENT_HELP,
        parents=[get_segment_parser()],
    )
    parser.set_defaults(func=cmd_run_segment)

    return parser


def cmd_run_segment(args):
    from openst.segmentation.segment import _run_segment

    _run_segment(args)


SEGMENT_MERGE_HELP = "Merge two segmentation masks into one"
def get_segment_merge_parser():
    parser = argparse.ArgumentParser(
        description=SEGMENT_MERGE_HELP,
        allow_abbrev=False,
        add_help=False,
    )
    parser.add_argument(
        "--h5-in",
        type=str,
        default="",
        required=True,
        help="""If set, masks are loaded from the Open-ST h5 object (key in --mask-in),
             and segmentation is saved there (to the key under --mask-out)""",
    )
    parser.add_argument(
        "--mask-in",
        type=str,
        required=True,
        nargs=2,
        help="Path to the input segmentation masks - two of them!",
    )
    parser.add_argument(
        "--mask-out",
        type=str,
        required=True,
        help="Path (file or h5) where the merged mask will be saved",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=512,
        help="When prediction of the mask runs in separate chunks, this is the chunk square size (in pixels)",
    )
    parser.add_argument(
        "--chunked",
        action="store_true",
        help="If set, segmentation is computed at non-overlapping chunks of size '--chunk-size'",
    )
    parser.add_argument(
        "--num-workers",
        type=int,
        help="Number of CPU workers when --chunked is specified",
        required=False,
        default=-1,
    )
    return parser


def setup_segment_merge_parser(parent_parser):
    parser = parent_parser.add_parser(
        "segment_merge",
        help=SEGMENT_MERGE_HELP,
        parents=[get_segment_merge_parser()],
    )
    parser.set_defaults(func=cmd_run_segment_merge)

    return parser


def cmd_run_segment_merge(args):
    from openst.segmentation.segment_merge import _run_segment_merge

    _run_segment_merge(args)


IMAGE_STITCH_HELP = "Stitch image fields of view into a single image"
def get_image_stitch_parser():
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        add_help=False,
        description=IMAGE_STITCH_HELP,
    )

    parser.add_argument(
        "--image-indir",
        type=str,
        help="path to collection of images",
        required=True,
    )

    parser.add_argument(
        "--image-out",
        type=str,
        help="path where to save the image (must be a filename)",
        required=True,
    )

    parser.add_argument(
        "--imagej-bin",
        type=str,
        help="path to the ImageJ/Fiji executable. Must have the Grid Collection plugin available!",
        required=True,
    )

    parser.add_argument(
        "--microscope",
        type=str,
        help="microscope model or imaging strategy that was used for imaging",
        choices=["keyence"],
        required=True,
    )

    parser.add_argument(
        "--no-run",
        default=False,
        action="store_true",
        help="If set, do not run ImageJ, but return the command line",
    )

    parser.add_argument(
        "--rerun",
        default=False,
        action="store_true",
        help="If set, runs stitching even when the output file exists",
    )

    parser.add_argument(
        "--metadata",
        type=str,
        default="",
        help="Path where the metadata will be stored. If not specified, metadata is not saved.",
    )

    parser.add_argument(
        "--join-zstack-regex",
        type=str,
        default="",
        help="""When non empty, this specifies how to find the Z location 
             from the individual filename and will create a z-stack from single images.
             Example regex: 'Image_([0-9]*)_Z([0-9]*)_CH1.tif'""",
    )

    return parser


def setup_image_stitch_parser(parent_parser):
    parser = parent_parser.add_parser(
        "image_stitch",
        help=IMAGE_STITCH_HELP,
        parents=[get_image_stitch_parser()],
    )
    parser.set_defaults(func=cmd_run_image_stitch)

    return parser


def cmd_run_image_stitch(args):
    from openst.preprocessing.image_stitch import _run_image_stitch

    _run_image_stitch(args)


SPATIAL_STITCH_HELP = "Stitch Open-ST h5 tile objects into a single Open-ST h5 object"
def get_spatial_stitch_parser():
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description=SPATIAL_STITCH_HELP,
        add_help=False,
    )

    parser.add_argument(
        "--tiles",
        type=str,
        nargs="+",
        help="Path to h5 file, one per tile - separated by space",
        required=True,
    )

    parser.add_argument(
        "--tile-coordinates",
        type=str,
        help="Path to the coordinate system file",
        required=True,
    )

    parser.add_argument(
        "--tile-id",
        type=str,
        nargs="+",
        help="""(Mandatory if --tile-id-regex is not specified)
              Per tile file specified in --tiles, each entry in --tile-id maps a tile file to the
              tile IDs under the first column of the --tile-coordinates file.""",
        default=None,
    )

    parser.add_argument(
        "--tile-id-regex",
        type=str,
        help="""(Mandatory if --tile-id is not specified)
              "Regex to find tile id in file names, instead of specifying a list in --tile-id""",
        default=DEFAULT_REGEX_TILE_ID,
    )

    parser.add_argument(
        "--h5-out",
        type=str,
        help="Where the stitched spatial object will be written to",
        required=True,
    )

    parser.add_argument(
        "--tile-id-key",
        type=str,
        help="""Key of the h5 file (under /obs) where tile IDs are stored.
                If != 'tile_id', a new categorical column of this name will be generated for consistency.""",
        default="tile_id",
    )

    parser.add_argument(
        "--merge-output",
        type=str,
        help='how to merge tiles, can be "same", "unique", "first", "only"',
        choices=["same", "unique", "first", "only"],
        default="same",
    )

    parser.add_argument(
        "--join-output",
        type=str,
        help='how to join tiles, can be "inner", "outer"',
        choices=["inner", "outer"],
        default="outer",
    )

    parser.add_argument(
        "--no-reset-index",
        default=False,
        action="store_true",
        help="""If set, do not reset the obs_name index of the combined spatial object
                as 'obs_name:<tile_id_key>'; keep original 'obs_name'""",
    )

    parser.add_argument(
        "--no-transform",
        default=False,
        action="store_true",
        help="If set, spatial coordinates are not transformed - just combine tiles into a single spatial object",
    )

    parser.add_argument(
        "--metadata",
        type=str,
        default="",
        help="(Optional) Path where the metadata will be stored. If not specified, metadata is not saved.",
    )

    return parser


def setup_spatial_stitch_parser(parent_parser):
    parser = parent_parser.add_parser(
        "spatial_stitch",
        help=SPATIAL_STITCH_HELP,
        parents=[get_spatial_stitch_parser()],
    )
    parser.set_defaults(func=cmd_run_spatial_stitch)

    return parser


def cmd_run_spatial_stitch(args):
    from openst.preprocessing.spatial_stitch import _run_spatial_stitch

    _run_spatial_stitch(args)

IMAGE_PREPROCESS_HELP = "Restoration of imaging data with CUT model (as in Open-ST paper)"
def get_image_preprocess_parser():
    parser = argparse.ArgumentParser(
        description=IMAGE_PREPROCESS_HELP,
        allow_abbrev=False,
        add_help=False,
    )
    parser.add_argument(
        "--h5-in",
        type=str,
        default="",
        help="""If set, image is loaded from the Open-ST h5 object (key in --image-in),
             and retored image is saved there (to the key --image-out)""",
    )
    parser.add_argument(
        "--image-in",
        type=str,
        default="uns/spatial/staining_image",
        help="Key or path to the input image",
    )
    parser.add_argument(
        "--image-out",
        type=str,
        default="uns/spatial/staining_image_restored",
        help="Key or path where the restored image will be written into",
    )
    parser.add_argument(
        "--tile-size-px",
        type=int,
        default=512,
        help="The input image is split into squared tiles of side `--tile-size-px`, for inference."+ 
        "Larger values avoid boundary effects, but require more memory.",
    )
    parser.add_argument("--model", type=str, default="HE_CUT_rajewsky", help="CUT model used for image restoration")
    parser.add_argument(
        "--device",
        type=str,
        default="cpu",
        choices=["cpu", "cuda"],
        help="Device used to run CUT restoration model. Can be ['cpu', 'cuda']",
    )
    parser.add_argument(
        "--num-workers",
        type=int,
        default=-1,
        help="Number of CPU workers for parallel processing",
    )

    return parser


def setup_image_preprocess_parser(parent_parser):
    parser = parent_parser.add_parser(
        "image_preprocess",
        help=IMAGE_PREPROCESS_HELP,
        parents=[get_image_preprocess_parser()],
    )
    parser.set_defaults(func=cmd_run_image_preprocess)

    return parser


def cmd_run_image_preprocess(args):
    from openst.preprocessing.image_preprocess import _run_image_preprocess

    _run_image_preprocess(args)

BARCODE_PREPROCESSING_HELP = "Convert spatial barcode raw data into tabular files with barcodes and spatial coordinates"
def get_barcode_preprocessing_parser():
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        add_help=False,
        description=BARCODE_PREPROCESSING_HELP,
    )
    parser.add_argument("--fastq-in", type=str, required=True, help="Path to the fastq file")
    parser.add_argument(
        "--tilecoords-out",
        type=str,
        required=True,
        help="Directory where output files will be written to",
    )
    parser.add_argument(
        "--out-suffix",
        type=str,
        required=True,
        help="Suffix added to the name of the output files (i.e., extension)",
    )
    parser.add_argument(
        "--out-prefix",
        type=str,
        default="",
        help="(Optional) Prefix added to the name of the output files",
    )
    parser.add_argument(
        "--crop-seq",
        type=str,
        default=":",
        help=")Optional) A 'python-style' slice, used to crop input sequences",
    )
    parser.add_argument(
        "--rev-comp",
        action="store_true",
        help="(Optional) Apply reverse complementary after sequence cropping",
    )
    parser.add_argument(
        "--single-tile",
        action="store_true",
        help="(Optional) set if it is guarranteed that the input .fastq(.gz) file contains only a tile",
    )
    parser.add_argument(
        "--unsorted",
        action="store_true",
        help="(Optional) set when file is unsorted respect to tiles; might be slower",
    )

    return parser


def setup_barcode_preprocessing_parser(parent_parser):
    parser = parent_parser.add_parser(
        "barcode_preprocessing",
        help=BARCODE_PREPROCESSING_HELP,
        parents=[get_barcode_preprocessing_parser()],
    )
    parser.set_defaults(func=cmd_run_barcode_preprocessing)

    return parser


def cmd_run_barcode_preprocessing(args):
    from openst.preprocessing.barcode_preprocessing import _run_barcode_preprocessing

    _run_barcode_preprocessing(args)

REPORT_HELP = "Generate HTML reports from metadata files (json)"
def get_report_parser():
    parser = argparse.ArgumentParser(
        description=REPORT_HELP,
        allow_abbrev=False,
        add_help=False,
    )

    parser.add_argument(
        "--metadata",
        type=str,
        required=True,
        help="Path to the metadata file (json)",
    )
    parser.add_argument(
        "--html-out",
        type=str,
        required=True,
        help="Path where the output HTML file will be created",
    )
    return parser


def setup_report_parser(parent_parser):
    parser = parent_parser.add_parser(
        "report",
        help=REPORT_HELP,
        parents=[get_report_parser()],
    )
    parser.set_defaults(func=cmd_run_report)

    return parser


def cmd_run_report(args):
    from openst.metadata.report import _run_report

    _run_report(args)

TRANSCRIPT_ASSIGN_HELP = "Aggregate transcripts into segmented cells"
def get_transcript_assign_parser():
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        add_help=False,
        description=TRANSCRIPT_ASSIGN_HELP,
    )

    parser.add_argument(
        "--h5-in",
        type=str,
        help="Path to an already aligned Open-ST h5 object",
        required=True,
    )

    parser.add_argument(
        "--mask-in",
        type=str,
        help="""Path to image mask - a key in the Open-ST h5 object.
                Or, can be an image stored separately in the filesystem (when --mask-from-file is specified)
                Image data and ST coordinates must be pairwise aligned
                (implicit for the case of RNA-based segmentation)""",
        required=True,
    )

    parser.add_argument(
        "--spatial-key",
        type=str,
        help="""Key in the Open-ST h5 object where the aligned coordinates are stored, 
        e.g. 'spatial_pairwise_aligned_coarse' (after using 'openst pairwise_aligner')""",
        required=True,
    )

    parser.add_argument(
        "--h5-out",
        type=str,
        help="Path where the segmented Open-ST h5 object will be written into",
        required=True,
    )

    parser.add_argument(
        "--mask-from-file",
        default=False,
        action="store_true",
        help="If set, the image mask is loaded from an external file",
    )

    parser.add_argument(
        "--max-image-pixels",
        type=int,
        default=933120000,
        help="Upper bound for number of pixels in the images (prevents exception when opening very large images)",
    )

    parser.add_argument(
        "--shuffle-umi",
        default=False,
        action="store_true",
        help="If set, UMI locations will be shuffled. This can be used as a baseline for feature selection.",
    )

    parser.add_argument(
        "--metadata",
        type=str,
        default="",
        help="""Path where the metadata will be stored.
        If not specified, metadata is not saved.
        Warning: a report (via openst report) cannot be generated without metadata!""",
    )

    return parser


def setup_transcript_assign_parser(parent_parser):
    parser = parent_parser.add_parser(
        "transcript_assign",
        help=TRANSCRIPT_ASSIGN_HELP,
        parents=[get_transcript_assign_parser()],
    )
    parser.set_defaults(func=cmd_run_transcript_assign)

    return parser


def cmd_run_transcript_assign(args):
    from openst.alignment.transcript_assign import _run_transcript_assign

    _run_transcript_assign(args)

APPLY_TRANSFORM_HELP = "Apply a precomputed transformation matrix to the specified coordinates of an Open-ST h5 object"
def get_apply_transform_parser():
    parser = argparse.ArgumentParser(
        description=APPLY_TRANSFORM_HELP,
        allow_abbrev=False,
        add_help=False,
    )
    parser.add_argument(
        "--keypoints-in",
        type=str,
        required=True,
        help="Path to the json file containing keypoints",
    )
    parser.add_argument(
        "--h5-in",
        type=str,
        required=True,
        help="Path to the input h5ad file containing spatial coordinates",
    )
    parser.add_argument(
        "--per-tile",
        action="store_true",
        help="""(Optional) If set, transformations are applied per tile, from their keypoints. 
                Otherwise, a single transform is computed for all tiles.""",
    )
    parser.add_argument(
        "--spatial-key-in",
        type=str,
        default="obsm/spatial_pairwise_aligned_coarse",
        help="""Key of the Open-ST h5 object where the input spatial coordinates are read from""",
    )
    parser.add_argument(
        "--spatial-key-out",
        type=str,
        default="obsm/spatial_pairwise_aligned_fine",
        help="""Key of the Open-ST h5 object where the transformed spatial coordinates are written into""",
    )
    return parser


def setup_apply_transform_parser(parent_parser):
    parser = parent_parser.add_parser(
        "apply_transform",
        help=APPLY_TRANSFORM_HELP,
        parents=[get_apply_transform_parser()],
    )
    parser.set_defaults(func=cmd_run_apply_transform)

    return parser


def cmd_run_apply_transform(args):
    from openst.alignment.apply_transform import _run_apply_transform

    _run_apply_transform(args)


MANUAL_PAIRWISE_ALIGNER_HELP = "GUI for manual alignment of Open-ST data"
def get_manual_pairwise_aligner_parser():
    parser = argparse.ArgumentParser(
        description=MANUAL_PAIRWISE_ALIGNER_HELP,
        allow_abbrev=False,
        add_help=False,
    )
    parser.add_argument(
        "--h5-in",
        type=str,
        default="",
        help="Path to the input h5ad file containing spatial coordinates",
    )
    parser.add_argument(
        "--spatial-key",
        type=str,
        default="",
        help="Path in the h5ad file to the spatial coordinates",
    )
    parser.add_argument(
        "--image-key",
        type=str,
        default="",
        help="Path in the h5ad file to the image",
    )
    return parser


def setup_manual_pairwise_aligner_parser(parent_parser):
    parser = parent_parser.add_parser(
        "manual_pairwise_aligner",
        help=MANUAL_PAIRWISE_ALIGNER_HELP,
        parents=[get_manual_pairwise_aligner_parser()],
    )
    parser.set_defaults(func=cmd_run_manual_pairwise_aligner)

    return parser


def cmd_run_manual_pairwise_aligner(args):
    from openst.alignment.manual_pairwise_aligner import _run_manual_pairwise_aligner

    _run_manual_pairwise_aligner(args)

PAIRWISE_ALIGNER_HELP = "Automatic pairwise alignment of transcript locations to imaging data"
def get_pairwise_aligner_parser():
    parser = argparse.ArgumentParser(
        description=PAIRWISE_ALIGNER_HELP,
        allow_abbrev=False,
        add_help=False,
    )

    required = parser.add_argument_group('Data (required)')
    required.add_argument(
        "--h5-in",
        type=str,
        required=True,
        help="Path to the merged Open-ST h5 object containing spatial coordinates and images",
    )

    optional_data = parser.add_argument_group('Data (optional)')
    optional_data.add_argument(
        "--image-in",
        type=str,
        default="uns/spatial/staining_image",
        help="Key to the image used as the 'destination' during pairwise alignment",
    )
    optional_data.add_argument(
        "--metadata",
        type=str,
        default="",
        help="Path where the metadata will be stored. If not specified, metadata is not saved.",
    )

    coarse_params = parser.add_argument_group('Coarse registration parameters')
    coarse_params.add_argument(
        "--only-coarse",
        action="store_true",
        help="If selected, only the coarse alignment stage will run",
    )
    coarse_params.add_argument(
        "--rescale-factor-coarse",
        type=int,
        default=20,
        help="Rescaling factor for the input image (1:factor), used during coarse pairwise alignment",
    )
    coarse_params.add_argument(
        "--threshold-counts-coarse",
        type=int,
        default=1,
        help="""Only spatial coordinates with counts larger than this number
        will be kept for pseudoimage rendering during coarse alignment""",
    )
    coarse_params.add_argument(
        "--pseudoimage-size-coarse",
        type=int,
        default=500,
        help="Size (in pixels) of the pseudoimage during coarse alignment.",
    )
    coarse_params.add_argument(
        "--ransac-coarse-min-samples",
        type=int,
        default=3,
        help="'min_samples' parameter of RANSAC, during coarse registration",
    )
    coarse_params.add_argument(
        "--ransac-coarse-residual-threshold",
        type=float,
        default=2,
        help="'residual_threshold' parameter of RANSAC, during coarse registration",
    )
    coarse_params.add_argument(
        "--ransac-coarse-max-trials",
        type=int,
        default=2,
        help="Times RANSAC will run (x1000 iterations) during coarse registration",
    )

    fine_params = parser.add_argument_group('Fine registration parameters')
    fine_params.add_argument(
        "--rescale-factor-fine",
        type=int,
        default=10,
        help="Rescaling factor for the input image (1:factor), used during fine pairwise alignment",
    )
    fine_params.add_argument(
        "--gaussian-sigma-fine",
        type=int,
        default=2,
        help="Gaussian blur used on all modalities during fine registration",
    )
    fine_params.add_argument(
        "--threshold-counts-fine",
        type=int,
        default=0,
        help="""Only spatial coordinates with counts larger than this number
        will be kept for pseudoimage rendering during fine alignment""",
    )
    fine_params.add_argument(
        "--pseudoimage-size-fine",
        type=int,
        default=2000,
        help="Size (in pixels) of the pseudoimage during fine alignment.",
    )
    
    fine_params.add_argument(
        "--ransac-fine-min-samples",
        type=int,
        default=3,
        help="'min_samples' parameter of RANSAC, during fine registration",
    )
    fine_params.add_argument(
        "--ransac-fine-residual-threshold",
        type=float,
        default=2,
        help="'residual_threshold' parameter of RANSAC, during fine registration",
    )
    fine_params.add_argument(
        "--ransac-fine-max-trials",
        type=int,
        default=1,
        help="Times RANSAC will run (x1000 iterations) during fine registration",
    )
    fine_params.add_argument(
        "--min-matches",
        type=int,
        default=50,
        help="Minimum number of matching keypoints between modalities during fine alignment",
    )

    image_preproc = parser.add_argument_group('Image preprocessing parameters')
    image_preproc.add_argument(
        "--mask-tissue",
        action="store_true",
        help="Tissue (imaging modality) is masked from the background for the feature detection",
    )
    fine_params.add_argument(
        "--mask-gaussian-sigma",
        type=int,
        default=5,
        help="The gaussian blur sigma used during the isolation of the tissue on the HE (preprocessing)",
    )
    image_preproc.add_argument(
        "--keep-black-background",
        action="store_true",
        help="Whether to set the background of the imaging modalities to white, after tissue masking",
    )
    
    model_params = parser.add_argument_group('Feature model parameters')
    model_params.add_argument(
        "--feature-matcher",
        type=str,
        default="LoFTR",
        choices=["LoFTR", "SIFT", "KeyNet"],
        help="Feature matching algorithm",
    )

    compu_params = parser.add_argument_group('Computational parameters')
    compu_params.add_argument(
        "--num-workers",
        type=int,
        default=1,
        help="Number of CPU workers for parallel processing",
    )
    compu_params.add_argument(
        "--device",
        type=str,
        default="cpu",
        choices=["cpu", "cuda"],
        help="Device used to run feature matching model. Can be ['cpu', 'cuda']",
    )
    return parser


def setup_pairwise_aligner_parser(parent_parser):
    parser = parent_parser.add_parser(
        "pairwise_aligner",
        help=PAIRWISE_ALIGNER_HELP,
        parents=[get_pairwise_aligner_parser()],
    )
    parser.set_defaults(func=cmd_run_pairwise_aligner)

    return parser


def cmd_run_pairwise_aligner(args):
    from openst.alignment.pairwise_aligner import _run_pairwise_aligner

    _run_pairwise_aligner(args)

FROM_3D_REGISTRATION_HELP = "Convert Open-ST h5 objects for 3D registration of serial sections using STIM"
def get_from_3d_registration_parser():
    parser = argparse.ArgumentParser(
        description=FROM_3D_REGISTRATION_HELP,
        allow_abbrev=False,
        add_help=False,
    )
    parser.add_argument(
        "--n5-dirs",
        type=str,
        nargs="+",
        required=True,
        help="Path to the n5 directories created by STIM for each of the samples. Expects more than one argument.",
    )
    parser.add_argument(
        "--h5-files",
        type=str,
        nargs="+",
        required=True,
        help="""Path to the h5ad files used to create the STIM data.
        Expects more than one argument, in the same order as --n5-dirs.""",
    )
    parser.add_argument(
        "--images",
        type=str,
        nargs="+",
        help="""Path to the file to load pairwise aligned images from.
        In this case, it should have the same length as the --n5-dirs and --h5-files arguments.
        If --images-in-adata is provided, it is the path inside the h5ad files where to load the images from.
        In this case, there only one argument is provided.""",
    )
    parser.add_argument(
        "--images-in-adata",
        default=False,
        action="store_true",
        help="""If set, the pairwise-aligned image is loaded from the adata,
        at the internal path specified by '--images'""",
    )
    parser.add_argument(
        "--h5-out",
        type=str,
        required=True,
        help="""Path to output h5ad file that will be created.
        When --separate-images is specified, it will not contain images.""",
    )
    parser.add_argument(
        "--image-out",
        type=str,
        help="""Path to output h5ad file that will be created.
        When --separate-images is specified, it will not contain images.""",
    )
    parser.add_argument(
        "--save-images-separately",
        default=False,
        action="store_true",
        help="When images are provided, these will be saved separately",
    )
    parser.add_argument(
        "--rescale",
        type=float,
        default=1,
        help="Rescales (multiples) the coordinates by --rescale units",
    )
    parser.add_argument(
        "--downsample-iamge",
        type=int,
        default=1,
        help="""How much to downwsample the resolution of the final aligned image volume
        (recommended for very large images). Must be >= 1""",
    )
    parser.add_argument(
        "--merge-output",
        type=str,
        help='how to merge tiles, can be "same", "unique", "first", "only"',
        choices=["same", "unique", "first", "only"],
        default="same",
    )
    parser.add_argument(
        "--join-output",
        type=str,
        help='how to join tiles, can be "inner", "outer"',
        choices=["inner", "outer"],
        default="outer",
    )
    parser.add_argument(
        "--max-image-pixels",
        type=int,
        default=933120000,
        help="Upper bound for number of pixels in the images (prevents exception when opening very large images)",
    )
    return parser


def setup_from_3d_registration_parser(parent_parser):
    parser = parent_parser.add_parser(
        "from_3d_registration",
        help=FROM_3D_REGISTRATION_HELP,
        parents=[get_from_3d_registration_parser()],
    )
    parser.set_defaults(func=cmd_run_from_3d_registration)

    return parser


def cmd_run_from_3d_registration(args):
    from openst.threed.from_3d_registration import _run_from_3d_registration

    _run_from_3d_registration(args)

GET_3D_REGISTRATION_HELP = """Convert STIM output back to a single (aligned) Open-ST h5 object.
                            If available, pairwise-aligned image data is transformed, too."""
def get_to_3d_registration_parser():
    parser = argparse.ArgumentParser(
        description=GET_3D_REGISTRATION_HELP,
        allow_abbrev=False,
        add_help=False,
    )
    parser.add_argument(
        "--in-adata",
        type=str,
        required=True,
        help="Path to the input adata file. Expects a file with raw counts.",
    )
    parser.add_argument(
        "--output-dir",
        type=int,
        default=-1,
        help="Path to the directory where the .genes and .locations csv files will be stored.",
    )
    parser.add_argument(
        "--filter-umi-max",
        type=int,
        default=-1,
        help="Preserve cells with at most < --filter-umi-max UMIs",
    )
    parser.add_argument(
        "--rescale",
        type=float,
        default=1,
        help="Rescales (multiples) the coordinates by --rescale units",
    )
    parser.add_argument(
        "--genes",
        type=str,
        nargs="+",
        default=None,
        help="Exports the expression level for the specified genes",
    )
    return parser


def setup_to_3d_registration_parser(parent_parser):
    parser = parent_parser.add_parser(
        "to_3d_registration",
        help=GET_3D_REGISTRATION_HELP,
        parents=[get_to_3d_registration_parser()],
    )
    parser.set_defaults(func=cmd_run_to_3d_registration)

    return parser


def cmd_run_to_3d_registration(args):
    from openst.threed.to_3d_registration import _run_to_3d_registration

    _run_to_3d_registration(args)

FROM_SPACEMAKE_HELP = """Run openst commands using spacemake file structure"""
def get_from_spacemake_parser():
    parser = argparse.ArgumentParser(
        description=FROM_SPACEMAKE_HELP,
        allow_abbrev=False,
        add_help=False,
    )
    parser.add_argument(
        "--project-id",
        type=str,
        required=True,
        help="From spacemake's project_df, this is the project_id string",
    )
    parser.add_argument(
        "--sample-id",
        type=str,
        required=True,
        help="From spacemake's project_df, this is the sample_id string",
    )
    parser.add_argument(
        "--run-mode",
        type=str,
        default="",
        help="When a sample has multiple run_mode(s), you must specify one",
    )
    return parser


def setup_from_spacemake_parser(parent_parser, from_spacemake_parser):
    parser = parent_parser.add_parser(
        "from_spacemake",
        help=FROM_SPACEMAKE_HELP,
        parents=[from_spacemake_parser],
    )
    parser.set_defaults(func=cmd_run_from_spacemake)

    return parser

def cmd_run_from_spacemake(parser, args, unknown_args):
    from openst.utils.from_spacemake import _run_from_spacemake

    _run_from_spacemake(parser, args, unknown_args)


def setup_from_spacemake_subparsers(from_spacemake_subparsers):
    setup_image_stitch_parser(from_spacemake_subparsers)
    setup_spatial_stitch_parser(from_spacemake_subparsers)

    setup_pairwise_aligner_parser(from_spacemake_subparsers)
    setup_apply_transform_parser(from_spacemake_subparsers)
    setup_manual_pairwise_aligner_parser(from_spacemake_subparsers)

    setup_segment_parser(from_spacemake_subparsers)
    setup_segment_merge_parser(from_spacemake_subparsers)
    setup_transcript_assign_parser(from_spacemake_subparsers)
    
    setup_merge_modalities_parser(from_spacemake_subparsers)
    setup_pseudoimage_parser(from_spacemake_subparsers)
    setup_preview_parser(from_spacemake_subparsers)

def cmdline_args():
    parent_parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description="Computational tools for Open-ST data",
    )

    description = """
    Flow cell preprocessing
        barcode_preprocessing
                        Convert spatial barcode raw data into tabular files with barcodes and spatial coordinates
    
    Sample preprocessing
        image_stitch        Stitch image fields of view into a single image
        image_preprocess    Restoration of imaging data with CUT model (as in Open-ST paper)
        spatial_stitch      Stitch Open-ST h5 tile objects into a single Open-ST h5 object
    
    Pairwise alignment
        merge_modalities    merge_modalities locations (as points) and images of Open-ST data
        pairwise_aligner    Automatic pairwise alignment of transcript locations to imaging data
        apply_transform     Apply a precomputed transformation matrix to the specified coordinates of an Open-ST h5 object
        manual_pairwise_aligner
                            GUI for manual alignment of Open-ST data
        
    Segmentation & DGE creation    
        segment             Image (or pseudoimage)-based segmentation with cellpose and (optional) radial extension
        segment_merge       Merge two segmentation masks into one
        transcript_assign   Aggregate transcripts into segmented cells
    
    3D reconstruction
        to_3d_registration  Convert STIM output back to a single (aligned) Open-ST h5 object. If available, pairwise-aligned image data is
                            transformed, too.
        from_3d_registration
                            Convert Open-ST h5 objects for 3D registration of serial sections using STIM
    
    Visualization and extra                        
        pseudoimage         Generate pseudoimages of Open-ST RNA data and visualize using napari
        preview             Preview locations (as points) and images of Open-ST data
        report              Generate HTML reports from metadata files (json)        
        from_spacemake      Run openst commands using spacemake file structure
    """

    parent_parser_subparsers = parent_parser.add_subparsers(title="commands", dest="subcommand")
    parent_parser.add_argument(
    '--version',
    action = 'store_true')

    # TODO: do this iteratively
    setup_barcode_preprocessing_parser(parent_parser_subparsers)
    
    setup_image_stitch_parser(parent_parser_subparsers)
    setup_image_preprocess_parser(parent_parser_subparsers)
    setup_spatial_stitch_parser(parent_parser_subparsers)

    setup_merge_modalities_parser(parent_parser_subparsers)
    setup_pairwise_aligner_parser(parent_parser_subparsers)
    setup_apply_transform_parser(parent_parser_subparsers)
    setup_manual_pairwise_aligner_parser(parent_parser_subparsers)

    setup_segment_parser(parent_parser_subparsers)
    setup_segment_merge_parser(parent_parser_subparsers)
    setup_transcript_assign_parser(parent_parser_subparsers)
    
    setup_to_3d_registration_parser(parent_parser_subparsers)
    setup_from_3d_registration_parser(parent_parser_subparsers)
    
    setup_pseudoimage_parser(parent_parser_subparsers)
    setup_preview_parser(parent_parser_subparsers)
    setup_report_parser(parent_parser_subparsers)

    # TODO: add subparser for from_spacemake
    from_spacemake_parser = get_from_spacemake_parser()
    setup_from_spacemake_parser(parent_parser_subparsers, from_spacemake_parser)

    known, unknown = parent_parser.parse_known_args()

    return parent_parser, from_spacemake_parser, known, unknown


def cmdline_main():
    import importlib.metadata
    parser, from_spacemake_parser, args, unknown_args = cmdline_args()

    if args.version and args.subcommand is None:
        print(importlib.metadata.version('openst'))
        return 0
    else:
        del args.version

    if "func" in args and args.subcommand != "from_spacemake":
        logging.info(f"openst {args.subcommand} - running with the following parameters:")
        logging.info(args.__dict__)
        args.func(args)
    elif args.subcommand == "from_spacemake":
        logging.info(f"openst {args.subcommand} - running from spacemake directory")
        logging.info(args.__dict__)
        args.func(from_spacemake_parser, args, unknown_args)
    else:
        parser.print_help()
        return 0


if __name__ == "__main__":
    cmdline_main()
