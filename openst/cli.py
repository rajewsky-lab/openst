import argparse

from openst.alignment.manual_pairwise_aligner import \
    setup_manual_pairwise_aligner_parser
from openst.alignment.manual_pairwise_aligner_gui import \
    setup_manual_pairwise_aligner_gui_parser
from openst.alignment.pairwise_aligner import setup_pairwise_aligner_parser
from openst.alignment.transcript_assign import setup_transcript_assign_parser
from openst.metadata.report import setup_report_parser
from openst.preprocessing.barcode_preprocessing import \
    setup_barcode_preprocessing_parser
from openst.preprocessing.image_preprocess import setup_image_preprocess_parser
from openst.preprocessing.image_stitch import setup_image_stitch_parser
from openst.preprocessing.spatial_stitch import setup_spatial_stitch_parser
from openst.segmentation.segment import setup_segment_parser
from openst.segmentation.segment_merge import setup_segment_merge_parser
from openst.threed.from_3d_registration import \
    setup_from_3d_registration_parser
from openst.threed.to_3d_registration import setup_to_3d_registration_parser


def cmdline_args():
    parent_parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description="openst: computational tools of Open-ST",
    )
    parent_parser_subparsers = parent_parser.add_subparsers(help="sub-command help", dest="subcommand")

    # TODO: do this iteratively
    # create the parser for the "pairwise_aligner" command
    setup_pairwise_aligner_parser(parent_parser_subparsers)
    # create the parser for the "manual_pairwise_aligner" command
    setup_manual_pairwise_aligner_parser(parent_parser_subparsers)
    # create the parser for the "manual_pairwise_aligner" command
    setup_manual_pairwise_aligner_gui_parser(parent_parser_subparsers)
    # create the parser for the "report" command
    setup_report_parser(parent_parser_subparsers)
    # create the parser for the "segment" command
    setup_segment_parser(parent_parser_subparsers)
    # create the parser for the "segment_merge" command
    setup_segment_merge_parser(parent_parser_subparsers)
    # create the parser for the "spatial_stitch" command
    setup_spatial_stitch_parser(parent_parser_subparsers)
    # create the parser for the "image_preprocess" command
    setup_image_preprocess_parser(parent_parser_subparsers)
    # create the parser for the "image_stitch" command
    setup_image_stitch_parser(parent_parser_subparsers)
    # create the parser for the "transcript_assign" command
    setup_transcript_assign_parser(parent_parser_subparsers)
    # create the parser for the "to_3d_registration" command
    setup_to_3d_registration_parser(parent_parser_subparsers)
    # create the parser for the "from_3d_registration" command
    setup_from_3d_registration_parser(parent_parser_subparsers)
    # create the parser for the "barcode_preprocessing" command
    setup_barcode_preprocessing_parser(parent_parser_subparsers)

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
