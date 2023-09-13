import os

import sys

class BaseMetadata:
    def __init__(self, args):
        self.metadata_type = None
        self.cmdline = " ".join(sys.argv)
        self.args = args

    def add_alignment_result(self, alignment_result):
        self.alignment_results.append(alignment_result)