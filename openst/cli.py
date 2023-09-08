import argparse

from openst.diffusion import setup_diffusion_parser


def cmdline_args():
    # create the top-level parser
    parent_parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description="Computational analysis of the open-ST pipeline",
    )
    parent_parser_subparsers = parent_parser.add_subparsers(help="sub-command help", dest="subcommand")

    # create the parser for the "diffusion" command
    setup_diffusion_parser(parent_parser_subparsers)
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
