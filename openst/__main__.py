import logging

from openst.cli import cmdline_main

logging.basicConfig(format="%(levelname)s %(asctime)s - %(message)s", level=logging.INFO)

if __name__ == "__main__":
    cmdline_main()
