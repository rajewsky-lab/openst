import argparse
import json
import logging
import os

from jinja2 import Template

from openst.utils.file import check_directory_exists, check_file_exists

absolute_path = os.path.dirname(__file__)


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
    parser.set_defaults(func=_run_report)

    return parser


def generate_html_report(json_data, template_file):
    with open(template_file, "r") as template_data:
        template_content = template_data.read()
        template = Template(template_content)

    html_report = template.render(
        program={"args": json_data["args"], "command_line": json_data["cmdline"]},
        metadata=json_data["metadata"]["alignment_results"],
    )

    return html_report


def _run_report(args):
    logging.info("openst report - running with the following parameters:")
    logging.info(args.__dict__)

    # Check input and output data
    check_file_exists(args.metadata)

    if not check_directory_exists(args.html_out):
        raise FileNotFoundError("Parent directory for --html-out does not exist")

    if check_file_exists(args.html_out, exception=False):
        logging.warn(f"The output file {args.html_out} exists; will be overwritten!")

    # paths provided by the user
    json_file = args.metadata
    output_file = args.html_out

    with open(json_file, "r") as json_data:
        json_data = json.load(json_data)

    # choose depending on the json-associated class

    template_file = os.path.join(absolute_path, f"templates/{json_data['metadata_type']}.html")

    html_report = generate_html_report(json_data, template_file)

    with open(output_file, "w") as output:
        output.write(html_report)


if __name__ == "__main__":
    args = get_report_parser().parse_args()
    _run_report(args)
