import logging
import os

from openst.cli import cmdline_main


def run_openst():
    OPENST_DEBUG = os.environ.get('OPENST_DEBUG', '0')
    _level = logging.DEBUG if OPENST_DEBUG == '1' else logging.INFO
    logging.basicConfig(format="%(levelname)s %(asctime)s - %(message)s", level=_level)
    cmdline_main()


if __name__ == "__main__":
    run_openst()
