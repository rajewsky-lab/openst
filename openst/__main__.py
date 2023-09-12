import logging

from openst.cli import cmdline_main


def run_openst():
    logging.basicConfig(format="%(levelname)s %(asctime)s - %(message)s", level=logging.INFO)
    cmdline_main()


if __name__ == "__main__":
    run_openst()
