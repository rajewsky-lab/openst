import argparse

from openst.threed.from_3d_registration import \
    setup_from_3d_registration_parser
from openst.threed.to_3d_registration import setup_to_3d_registration_parser

DEFAULT_REGEX_TILE_ID = "(L[1-4][a-b]_tile_[1-2][0-7][0-9][0-9])"

def get_pseudoimage_parser():
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="generate and visualize pseudoimage of Open-ST RNA data",
        allow_abbrev=False,
        add_help=False,
    )
    # Input
    parser.add_argument(
        "--adata",
        type=str,
        required=True,
        help="Necessary to create the pseudoimage",
    )

    # RNA density-based pseudoimage
    parser.add_argument(
        "--spatial-coord-key",
        type=str,
        default='obsm/spatial',
        help="Path to the spatial coordinates inside the AnnData object (e.g., 'obsm/spatial')"
    )
    parser.add_argument(
        "--input-resolution",
        type=float,
        default=1,
        help="""Spatial resolution of the input coordinates (retrieved from --spatial-coord-key).
              If it is in microns, leave as 1. If it is in pixels, specify the pixel to micron conversion factor."""
    )
    parser.add_argument(
        "--render-scale",
        type=float,
        default=2,
        help="Size of bins for computing the binning (in microns). For Open-ST v1, we recommend a value of 2."
    )
    parser.add_argument(
        "--render-sigma",
        type=float,
        default=1,
        help="Smoothing factor applied to the RNA pseudoimage (higher values lead to smoother images)"
    )
    parser.add_argument(
        "--output-resolution",
        type=float,
        default=0.6,
        help="Final resolution (micron/pixel) for the segmentation mask."
    )
    return parser

def setup_pseudoimage_parser(parent_parser):
    """setup_pseudoimage_parser"""
    parser = parent_parser.add_parser(
        "pseudoimage",
        help="generate and visualize pseudoimage of Open-ST RNA data",
        parents=[get_pseudoimage_parser()],
    )
    parser.set_defaults(func=cmd_run_pseudoimage_visualizer)

    return parser

def cmd_run_pseudoimage_visualizer(args):
    from openst.utils.pseudoimage import _run_pseudoimage_visualizer
    _run_pseudoimage_visualizer(args)


def get_preview_parser():
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="preview pseudoimage and images of Open-ST data",
        allow_abbrev=False,
        add_help=False,
    )
    # Input
    parser.add_argument(
        "--adata",
        type=str,
        required=True,
        help="Necessary to create the pseudoimage",
    )

    # RNA density-based pseudoimage
    parser.add_argument(
        "--spatial-coord-keys",
        type=str,
        nargs="+",
        default=None,
        help="""Path to the spatial coordinates inside the AnnData object (e.g., 'obsm/spatial').
                Can be one or many (separated by space)"""
    )

    # Staining image
    parser.add_argument(
        "--image-keys",
        type=str,
        nargs="+",
        default=None,
        help="""Path to the image to be visualized.
              Can be one or many (separated by space)"""
    )

    # Resampling before previsualizing
    parser.add_argument(
        "--spatial-coord-resampling",
        type=int,
        nargs="+",
        default=[1],
        help="""Will load every n-th point. Can be one (same for all spatial-coords)
                or many (1-to-1 mapping to the spatial-coord list)"""
    )
    parser.add_argument(
        "--image-resampling",
        type=int,
        nargs="+",
        default=[1],
        help="""Will load every n-th pixel. Can be one (same for all images)
                or many (1-to-1 mapping to the image list)"""
    )

    return parser

def setup_preview_parser(parent_parser):
    """setup_preview_parser"""
    parser = parent_parser.add_parser(
        "preview",
        help="prepreview the dataset with napari",
        parents=[get_preview_parser()],
    )
    parser.set_defaults(func=cmd_run_preview)

    return parser

def cmd_run_preview(args):
    from openst.utils.preview import _run_preview
    _run_preview(args)

def get_segment_parser():
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="segmentation of Open-ST imaging data with cellpose",
        allow_abbrev=False,
        add_help=False,
    )

    # Data
    parser.add_argument(
        "--image-in",
        type=str,
        help="Path to the input image.",
    )
    parser.add_argument(
        "--rna-segment",
        action="store_true",
        help="""Performs segmentation based on local RNA density pseudoimages from sequencing data,
              instead of using a staining image. 
              This assumes coordinates in microns (can be transformed with --rna-segment-input-resolution)"""
    )
    parser.add_argument(
        "--output-mask",
        type=str,
        required=True,
        help="Path to the output file where the mask will be stored",
    )
    parser.add_argument(
        "--adata",
        type=str,
        default="",
        help="When specified, staining image is loaded from adata (from --image-in), and segmentation is saved there (to --output-mask)",
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
        help="When specified, segmentation is computed at non-overlapping chunks of size '--chunk-size'",
    )
    parser.add_argument(
        "--max-image-pixels",
        type=int,
        default=933120000,
        help="Upper bound for number of pixels in the images (prevents exception when opening very large images)",
    )
    parser.add_argument(
        "--gpu",
        action="store_true",
        help="When specified, a GPU will be used for prediction",
    )

    # Postprocessing
    parser.add_argument(
        "--dilate-px",
        type=int,
        help="Pixels the outlines of the segmentation mask will be extended",
        required=False,
        default=0,
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
        default='obsm/spatial',
        help="Path to the spatial coordinates inside the AnnData object (e.g., 'obsm/spatial')"
    )
    parser.add_argument(
        "--rna-segment-input-resolution",
        type=float,
        default=1,
        help="""Spatial resolution of the input coordinates (retrieved from --rna-segment-spatial-coord-key).
              If it is in microns, leave as 1. If it is in pixels, specify the pixel to micron conversion factor."""
    )
    parser.add_argument(
        "--rna-segment-render-scale",
        type=float,
        default=2,
        help="Size of bins for computing the binning (in microns). For Open-ST v1, we recommend a value of 2."
    )
    parser.add_argument(
        "--rna-segment-render-sigma",
        type=float,
        default=1,
        help="Smoothing factor applied to the RNA pseudoimage (higher values lead to smoother images)"
    )
    parser.add_argument(
        "--rna-segment-output-resolution",
        type=float,
        default=0.6,
        help="Final resolution (micron/pixel) for the segmentation mask."
    )
    

    # General
    parser.add_argument(
        "--num-workers",
        type=int,
        help="Number of parallel workers when --chunked is specified",
        required=False,
        default=-1,
    )
    parser.add_argument(
        "--metadata-out",
        type=str,
        default="",
        help="""Path where the metadata will be stored.
        If not specified, metadata is not saved.
        Warning: a report (via openst report) cannot be generated without metadata!""",
    )
    return parser


def setup_segment_parser(parent_parser):
    """setup_segment_parser"""
    parser = parent_parser.add_parser(
        "segment",
        help="segmentation of Open-ST imaging data with cellpose",
        parents=[get_segment_parser()],
    )
    parser.set_defaults(func=cmd_run_segment)

    return parser


def cmd_run_segment(args):
    from openst.segmentation.segment import _run_segment
    _run_segment(args)


def get_segment_merge_parser():
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="segmentation of Open-ST imaging data with cellpose",
        allow_abbrev=False,
        add_help=False,
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
        "--adata",
        type=str,
        default='',
        help="When specified, masks are loaded from adata (--mask-in), and segmentation is saved there (to --mask-out)",
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
        help="When specified, segmentation is computed at non-overlapping chunks of size '--chunk-size'",
    )
    parser.add_argument(
        "--max-image-pixels",
        type=int,
        default=933120000,
        help="Upper bound for number of pixels in the images (prevents exception when opening very large images)",
    )
    parser.add_argument(
        "--num-workers",
        type=int,
        help="Number of parallel workers when --chunked is specified",
        required=False,
        default=-1,
    )
    return parser


def setup_segment_merge_parser(parent_parser):
    """setup_segment_merge_parser"""
    parser = parent_parser.add_parser(
        "segment_merge",
        help="merge segmentations from Open-ST data",
        parents=[get_segment_merge_parser()],
    )
    parser.set_defaults(func=cmd_run_segment_merge)

    return parser


def cmd_run_segment_merge(args):
    from openst.segmentation.segment_merge import _run_segment_merge
    _run_segment_merge(args)


def get_image_stitch_parser():
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        add_help=False,
        description="stitching various image tiles into a single image (wrapper for ImageJ)",
    )

    parser.add_argument(
        "--input-dir",
        type=str,
        help="path to collection of images",
        required=True,
    )

    parser.add_argument(
        "--output-image",
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
        help="When specified, do not run ImageJ, but return the command line",
    )

    parser.add_argument(
        "--metadata-out",
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
             Example regex: 'Image_([0-9]*)_Z([0-9]*)_CH1.tif'"""
    )

    return parser


def setup_image_stitch_parser(parent_parser):
    """setup_image_stitch_parser"""
    parser = parent_parser.add_parser(
        "image_stitch",
        help="stitching image tiles into a single image (wrapper for ImageJ)",
        parents=[get_image_stitch_parser()],
    )
    parser.set_defaults(func=cmd_run_image_stitch)

    return parser


def cmd_run_image_stitch(args):
    from openst.preprocessing.image_stitch import _run_image_stitch
    _run_image_stitch(args)


def get_spatial_stitch_parser():
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description="stitching spatial transcriptomics tiles into a common global coordinate system",
        add_help=False,
    )

    parser.add_argument(
        "--tiles",
        type=str,
        nargs="+",
        help="path to spatial.h5ad AnnData file, one per tile",
        required=True,
    )

    parser.add_argument(
        "--tile-id",
        type=str,
        nargs="+",
        help="list of tile id for the input files, same order as tiles."
        + "Must be specified when the filenames do not contain a tile id that can be parsed with --tile-id-regex",
        default=None,
    )

    parser.add_argument(
        "--tile-coordinates",
        type=str,
        help="name of tile collection",
        required=True,
    )

    parser.add_argument(
        "--output",
        type=str,
        help="path to output.h5ad AnnData file",
        required=True,
    )

    parser.add_argument(
        "--tile-id-regex",
        type=str,
        help="regex to find tile id in file names",
        default=DEFAULT_REGEX_TILE_ID,
    )

    parser.add_argument(
        "--tile-id-key",
        type=str,
        help="name of .obs variable where tile id are (will be) stored",
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
        help="""do not reset the obs_name index of the AnnData object
        as 'obs_name:<tile_id_key>'; keep original 'obs_name'""",
    )

    parser.add_argument(
        "--no-transform",
        default=False,
        action="store_true",
        help="do not transform the spatial coordinates of the AnnData object",
    )

    parser.add_argument(
        "--metadata-out",
        type=str,
        default="",
        help="Path where the metadata will be stored. If not specified, metadata is not saved.",
    )

    return parser


def setup_spatial_stitch_parser(parent_parser):
    """setup_spatial_stitch_parser"""
    parser = parent_parser.add_parser(
        "spatial_stitch",
        help="stitching STS tiles into a global coordinate system",
        parents=[get_spatial_stitch_parser()],
    )
    parser.set_defaults(func=cmd_run_spatial_stitch)

    return parser


def cmd_run_spatial_stitch(args):
    from openst.preprocessing.spatial_stitch import _run_spatial_stitch
    _run_spatial_stitch(args)


def get_image_preprocess_parser():
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="preprocess imaging data with CUT model (as in Open-ST paper)",
        allow_abbrev=False,
        add_help=False,
    )
    parser.add_argument('--input_img', type=str, required=True,
                        help='path to input image')
    parser.add_argument('--cut_dir', type=str, required=True,
                        help='path to CUT directory (to save patched images)')
    parser.add_argument('--tile_size_px', type=int, required=True,
                        help='size of the tile in pixels')
    parser.add_argument(
        "--device",
        type=str,
        default="cpu",
        choices=["cpu", "cuda"],
        help="Device used to run feature matching model. Can be ['cpu', 'cuda']",
    )
    parser.add_argument('--checkpoints_dir', type=str, default='./checkpoints', help='models are saved here')

    return parser


def setup_image_preprocess_parser(parent_parser):
    """setup_image_preprocess_parser"""
    parser = parent_parser.add_parser(
        "image_preprocess",
        help="convert fastq files into spatial barcode files (sequence and coordinates)",
        parents=[get_image_preprocess_parser()],
    )
    parser.set_defaults(func=cmd_run_image_preprocess)

    return parser

def cmd_run_image_preprocess(args):
    from openst.preprocessing.image_preprocess import _run_image_preprocess
    _run_image_preprocess(args)


def get_barcode_preprocessing_parser():
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        add_help=False,
        description="convert fastq files into spatial barcode files (sequence and coordinates)",
    )
    parser.add_argument("--in-fastq", type=str, required=True, help="path to the fastq file")
    parser.add_argument(
        "--out-path",
        type=str,
        required=True,
        help="folder where the output files will be generated",
    )
    parser.add_argument(
        "--out-suffix",
        type=str,
        required=True,
        help="where to write the output file. it works as the suffix when multiple tiles are generated",
    )
    parser.add_argument(
        "--out-prefix",
        type=str,
        default="",
        help="where to write the output file. it works as the prefix when multiple tiles are generated",
    )
    parser.add_argument(
        "--crop-seq",
        type=str,
        default=":",
        help="crop the input sequence, should be a valid python slice",
    )
    parser.add_argument(
        "--rev-comp",
        action="store_true",
        help="applies reverse complementary after sequence cropping",
    )
    parser.add_argument(
        "--single-tile",
        action="store_true",
        help="it is guarranteed that the input .fastq(.gz) file contains only a tile. Throw an error otherwise",
    )
    parser.add_argument(
        "--unsorted",
        action="store_true",
        help="supports that the file is unsorted respect to tiles. might be slower",
    )

    return parser


def setup_barcode_preprocessing_parser(parent_parser):
    """setup_barcode_preprocessing_parser"""
    parser = parent_parser.add_parser(
        "barcode_preprocessing",
        help="convert fastq files into spatial barcode files (sequence and coordinates)",
        parents=[get_barcode_preprocessing_parser()],
    )
    parser.set_defaults(func=cmd_run_barcode_preprocessing)

    return parser

def cmd_run_barcode_preprocessing(args):
    from openst.preprocessing.barcode_preprocessing import _run_barcode_preprocessing
    _run_barcode_preprocessing(args)


def get_report_parser():
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="openst report HTML generator from metadata files (json)",
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
    """setup_report_parser"""
    parser = parent_parser.add_parser(
        "report",
        help="openst report HTML generator from metadata files (json)",
        parents=[get_report_parser()],
    )
    parser.set_defaults(func=cmd_run_report)

    return parser

def cmd_run_report(args):
    from openst.metadata.report import _run_report
    _run_report(args)


def get_transcript_assign_parser():
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        add_help=False,
        description="openst transfer of transcripts to single cells using a pairwise-aligned segmentation mask",
    )

    parser.add_argument(
        "--adata",
        type=str,
        help="path to previously aligned spatial.h5ad AnnData file",
        required=True,
    )

    parser.add_argument(
        "--mask-in-adata",
        default=False,
        action="store_true",
        help="When specified, the image mask is loaded from the adata, at the internal path specified by '--mask'",
    )

    parser.add_argument(
        "--shuffle-umi",
        default=False,
        action="store_true",
        help="When specified, UMI locations will be shuffled. This can be used as a baseline for feature selection.",
    )

    parser.add_argument(
        "--mask",
        type=str,
        help="""path to image mask; must be in same coordinates as the obsm['spatial'] in the AnnData file.
        If --mask-in-adata, it is a path within the h5ad file""",
        required=True,
    )

    parser.add_argument(
        "--spatial-key",
        type=str,
        help="""obsm dataset for the aligned coordinates, e.g. 'spatial_pairwise_aligned_coarse'  or
        'spatial_pairwise_aligned_fine' if the data has been aligned with openst pairwise_aligner""",
        required=True,
    )

    parser.add_argument(
        "--output",
        type=str,
        help="path and filename for output file that will be generated",
        required=True,
    )

    parser.add_argument(
        "--max-image-pixels",
        type=int,
        default=933120000,
        help="Upper bound for number of pixels in the images (prevents exception when opening very large images)",
    )
    parser.add_argument(
        "--metadata-out",
        type=str,
        default="",
        help="""Path where the metadata will be stored.
        If not specified, metadata is not saved.
        Warning: a report (via openst report) cannot be generated without metadata!""",
    )

    return parser


def setup_transcript_assign_parser(parent_parser):
    """setup_transcript_assign_parser"""
    parser = parent_parser.add_parser(
        "transcript_assign",
        help="assign transcripts into previously aligned segmentation mask",
        parents=[get_transcript_assign_parser()],
    )
    parser.set_defaults(func=cmd_run_transcript_assign)

    return parser

def cmd_run_transcript_assign(args):
    from openst.alignment.transcript_assign import _run_transcript_assign
    _run_transcript_assign(args)


def get_manual_pairwise_aligner_parser():
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="openst pairwise alignment of two-dimensional spatial transcriptomics and imaging data",
        allow_abbrev=False,
        add_help=False,
    )
    parser.add_argument(
        "--keypoints-json",
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
        "--h5-out",
        type=str,
        default="",
        help="""Path where the h5ad file will be saved after alignment.
        If not indicated, data is written in place at --h5-in""",
    )
    parser.add_argument(
        "--per-tile",
        action="store_true",
        help="If selected, individual transformations per tile are estimated from they keypoints",
    )
    parser.add_argument(
        "--spatial-key-in",
        type=str,
        default="spatial_pairwise_aligned_coarse",
        help="""The name of the `obsm` variable where the transformed coordinates will be read from""",
    )
    parser.add_argument(
        "--spatial-key-out",
        type=str,
        default="spatial_pairwise_aligned_fine",
        help="""The name of the `obsm` variable where the transformed coordinates will be written""",
    )
    return parser


def setup_manual_pairwise_aligner_parser(parent_parser):
    """setup_manual_pairwise_aligner_parser"""
    parser = parent_parser.add_parser(
        "manual_pairwise_aligner",
        help="openst manual pairwise alignment of spatial transcriptomics and imaging data",
        parents=[get_manual_pairwise_aligner_parser()],
    )
    parser.set_defaults(func=cmd_run_manual_pairwise_aligner)

    return parser

def cmd_run_manual_pairwise_aligner(args):
    from openst.alignment.manual_pairwise_aligner import _run_manual_pairwise_aligner
    _run_manual_pairwise_aligner(args)


def setup_manual_pairwise_aligner_gui_parser(parent_parser):
    """setup_manual_pairwise_aligner_gui_parser"""
    parser = parent_parser.add_parser(
        "manual_pairwise_aligner_gui",
        help="GUI for openst manual pairwise alignment of spatial transcriptomics and imaging data",
    )
    parser.set_defaults(func=cmd_run_manual_pairwise_aligner_gui)

    return parser

def cmd_run_manual_pairwise_aligner_gui(args):
    from openst.alignment.manual_pairwise_aligner_gui import _run_manual_pairwise_aligner_gui
    _run_manual_pairwise_aligner_gui(args)


def get_pairwise_aligner_parser():
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="openst pairwise alignment of two-dimensional spatial transcriptomics and imaging data",
        allow_abbrev=False,
        add_help=False,
    )

    parser.add_argument(
        "--image-in",
        type=str,
        required=True,
        help="Path to the input image. This is treated as the 'destination' image during pairwise alignment",
    )
    parser.add_argument(
        "--h5-in",
        type=str,
        required=True,
        help="Path to the input h5ad file containing spatial coordinates",
    )
    parser.add_argument(
        "--metadata-out",
        type=str,
        default="",
        help="Path where the metadata will be stored. If not specified, metadata is not saved.",
    )
    parser.add_argument(
        "--h5-out",
        type=str,
        required=True,
        help="Path where the h5ad file will be saved after alignment",
    )
    parser.add_argument(
        "--save-image-in-h5",
        action="store_true",
        help="Whether the input image will be saved into the output h5ad file",
    )
    parser.add_argument(
        "--mask-tissue",
        action="store_true",
        help="Tissue (imaging modality) is masked from the background for the feature detection",
    )
    parser.add_argument(
        "--only-coarse",
        action="store_true",
        help="If selected, only the coarse alignment stage will run",
    )
    parser.add_argument(
        "--rescale-factor-coarse",
        type=int,
        default=20,
        help="Rescaling factor for the input image (1:factor), used during coarse pairwise alignment",
    )
    parser.add_argument(
        "--rescale-factor-fine",
        type=int,
        default=5,
        help="Rescaling factor for the input image (1:factor), used during fine pairwise alignment",
    )
    parser.add_argument(
        "--tissue-masking-gaussian-sigma",
        type=int,
        default=5,
        help="The gaussian blur sigma used during the isolation of the tissue on the HE (preprocessing)",
    )
    parser.add_argument(
        "--fine-registration-gaussian-sigma",
        type=int,
        default=2,
        help="Gaussian blur used on all modalities during fine registration",
    )
    parser.add_argument(
        "--keep-black-background",
        action="store_true",
        help="Whether to set the background of the imaging modalities to white, after tissue masking",
    )
    parser.add_argument(
        "--threshold-counts-coarse",
        type=int,
        default=1,
        help="""Only spatial coordinates with counts larger than this number
        will be kept for pseudoimage rendering during coarse alignment""",
    )
    parser.add_argument(
        "--threshold-counts-fine",
        type=int,
        default=0,
        help="""Only spatial coordinates with counts larger than this number
        will be kept for pseudoimage rendering during fine alignment""",
    )
    parser.add_argument(
        "--pseudoimage-size-coarse",
        type=int,
        default=4000,
        help="Size (in pixels) of the pseudoimage during coarse alignment.",
    )
    parser.add_argument(
        "--pseudoimage-size-fine",
        type=int,
        default=6000,
        help="Size (in pixels) of the pseudoimage during fine alignment.",
    )
    parser.add_argument(
        "--ransac-coarse-min-samples",
        type=int,
        default=3,
        help="'min_samples' parameter of RANSAC, during coarse registration",
    )
    parser.add_argument(
        "--ransac-coarse-residual-threshold",
        type=float,
        default=2,
        help="'residual_threshold' parameter of RANSAC, during coarse registration",
    )
    parser.add_argument(
        "--ransac-coarse-max-trials",
        type=int,
        default=50,
        help="Times RANSAC will run (x1000 iterations) during coarse registration",
    )
    parser.add_argument(
        "--ransac-fine-min-samples",
        type=int,
        default=3,
        help="'min_samples' parameter of RANSAC, during fine registration",
    )
    parser.add_argument(
        "--ransac-fine-residual-threshold",
        type=float,
        default=2,
        help="'residual_threshold' parameter of RANSAC, during fine registration",
    )
    parser.add_argument(
        "--ransac-fine-max-trials",
        type=int,
        default=50,
        help="Times RANSAC will run (x1000 iterations) during fine registration",
    )
    parser.add_argument(
        "--max-image-pixels",
        type=int,
        default=933120000,
        help="Upper bound for number of pixels in the images (prevents exception when opening very large images)",
    )
    parser.add_argument(
        "--n-threads",
        type=int,
        default=1,
        help="Number of CPU threads for parallel processing",
    )
    parser.add_argument(
        "--feature-matcher",
        type=str,
        default="LoFTR",
        choices=["LoFTR", "SIFT", 'KeyNet'],
        help="Feature matching algorithm",
    )
    parser.add_argument(
        "--fine-min-matches",
        type=int,
        default=50,
        help="Minimum number of matching keypoints between modalities during fine alignment",
    )
    parser.add_argument(
        "--fiducial-model",
        type=str,
        default="",
        help="Path to a object detection model (YOLO) to detect fiducial markers",
    )
    parser.add_argument(
        "--genes-coarse",
        nargs="+",
        type=str,
        default=None,
        help="Genes used for plotting the pseudoimage during the coarse alignment phase.",
    )
    parser.add_argument(
        "--genes-fine",
        nargs="+",
        type=str,
        default=None,
        help="Genes used for plotting the pseudoimage during the fine alignment phase.",
    )
    parser.add_argument(
        "--device",
        type=str,
        default="cpu",
        choices=["cpu", "cuda"],
        help="Device used to run feature matching model. Can be ['cpu', 'cuda']",
    )
    return parser


def setup_pairwise_aligner_parser(parent_parser):
    """setup_pairwise_aligner_parser"""
    parser = parent_parser.add_parser(
        "pairwise_aligner",
        help="openst pairwise alignment of two-dimensional spatial transcriptomics and imaging data",
        parents=[get_pairwise_aligner_parser()],
    )
    parser.set_defaults(func=cmd_run_pairwise_aligner)

    return parser

def cmd_run_pairwise_aligner(args):
    from openst.alignment.pairwise_aligner import _run_pairwise_aligner
    _run_pairwise_aligner(args)


def cmdline_args():
    parent_parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description="openst: computational tools of Open-ST",
    )
    parent_parser_subparsers = parent_parser.add_subparsers(help="sub-command help", dest="subcommand")

    # TODO: do this iteratively
    setup_pairwise_aligner_parser(parent_parser_subparsers)
    setup_manual_pairwise_aligner_parser(parent_parser_subparsers)
    setup_manual_pairwise_aligner_gui_parser(parent_parser_subparsers)
    setup_report_parser(parent_parser_subparsers)
    setup_segment_parser(parent_parser_subparsers)
    setup_segment_merge_parser(parent_parser_subparsers)
    setup_spatial_stitch_parser(parent_parser_subparsers)
    setup_image_preprocess_parser(parent_parser_subparsers)
    setup_image_stitch_parser(parent_parser_subparsers)
    setup_transcript_assign_parser(parent_parser_subparsers)
    setup_to_3d_registration_parser(parent_parser_subparsers)
    setup_from_3d_registration_parser(parent_parser_subparsers)
    setup_barcode_preprocessing_parser(parent_parser_subparsers)
    setup_pseudoimage_parser(parent_parser_subparsers)
    setup_preview_parser(parent_parser_subparsers)

    return parent_parser, parent_parser.parse_args()


def cmdline_main():
    parser, args = cmdline_args()

    if "func" in args:
        args.func(args)
    else:
        parser.print_help()
        return 0


if __name__ == "__main__":
    cmdline_main()
