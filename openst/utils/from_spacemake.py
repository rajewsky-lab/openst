import os
import pandas as pd
import yaml
import logging
import ast

from openst.utils.file import check_directory_exists, check_file_exists

config_path = "config.yaml"
project_df = "project_df.csv"
COMPATIBLE_SPACEMAKE_VERSIONS = [
    "0.5.2",
    "0.7.3",
    "0.7.5"
]
VALID_PROJECT_DF_COLUMNS = [
    "project_id",
    "sample_id",
    "puck_barcode_file_id",
    "sample_sheet",
    "species",
    "demux_barcode_mismatch",
    "demux_dir",
    "basecalls_dir",
    "R1",
    "R2",
    "longreads",
    "longread_signature",
    "investigator",
    "sequencing_date",
    "experiment",
    "puck_barcode_file",
    "run_mode",
    "barcode_flavor",
    "is_merged",
    "merged_from",
    "puck",
    "dge",
    "map_strategy"
]

SMK_DGE_TYPES = {
    ".exon",
    ".intron",
    ".all",
    ".Reads_exon",
    ".Reads_intron",
    ".Reads_all",
    "",
}

SMK_PROJECT = 'projects/{project_id}'
SMK_SAMPLE_PROCESSED = 'processed_data/{sample_id}'
SMK_SAMPLE_RAW = 'raw_data/{sample_id}'

SMK_DGE = os.path.join(SMK_PROJECT, SMK_SAMPLE_PROCESSED, 'illumina/complete_data/dge')
SMK_DGE_FILE = os.path.join(SMK_DGE, 'dge{dge_type}{dge_cleaned}{polyA_adapter_trimmed}{mm_included}.spatial_beads{tile}.h5ad')

SMK_IMAGES_PROCESSED = os.path.join(SMK_PROJECT, SMK_SAMPLE_PROCESSED, 'images')
SMK_IMAGES_PROCESSED_STITCHED = os.path.join(SMK_IMAGES_PROCESSED, 'Image_Stitched_Composite.tif')
SMK_IMAGES_RAW = os.path.join(SMK_PROJECT, SMK_SAMPLE_RAW, 'images')

SMK_MULTIMODAL_PROCESSED = os.path.join(SMK_PROJECT, SMK_SAMPLE_PROCESSED, 'multimodal')
SMK_MULTIMODAL_PROCESSED_DGE = os.path.join(SMK_MULTIMODAL_PROCESSED, 'stitched_spots.h5ad')

class SpacemakeError(Exception):
    def __init__(self, msg=None):
        self.msg = msg

    def __str__(self):
        msg = 'ERROR: ' + str(self.__class__.__name__) + '\n'

        if hasattr(self, 'msg') and self.msg is not None:
            msg += self.msg

        return msg

def get_sample_metadata(pdf, project_id, sample_id):
    _sample = pdf[(pdf['project_id'] == project_id) & (pdf['sample_id'] == sample_id)]
    if len(_sample) > 1:
        raise SpacemakeError(f"Entry project_id='{project_id}' and sample_id='{sample_id}' is duplicated")
    
    _sample = _sample.iloc[0].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) and "[" in x else x)
    return _sample.to_dict()
    

def is_project_df_valid(pdf):
    if len(set(pdf.columns).intersection(set(VALID_PROJECT_DF_COLUMNS))) != len(VALID_PROJECT_DF_COLUMNS):
        return False
    
    return True

def is_config_valid(config):
    return True
    
def load_project_df():
    pdf = pd.read_csv(project_df)
    if not is_project_df_valid(pdf):
        raise SpacemakeError("The current 'project_df.csv' is not valid")
    
    return pdf

def load_config():
    with open('config.yaml') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    is_config_valid(config)

    return config

def get_config_from_sample(sample_metadata, config, run_mode=None):
    sample_config = {}

    # Puck
    sample_config['puck'] = config['pucks'][sample_metadata['puck']]

    # Run mode options
    if run_mode is None and len(sample_metadata['run_mode']) > 1:
        raise SpacemakeError(f"""There are various 'run_modes' for this sample:
                             \t{sample_metadata['run_mode']}
                              Specify only one using the '--run-mode' argument""")
    elif run_mode is None:
        run_mode = sample_metadata['run_mode'][0]

    sample_config['run_mode'] = {}
    sample_config['run_mode']['dge_type'] = ".all" if config['run_modes'][run_mode]['count_intronic_reads'] else ".exon"
    sample_config['run_mode']['dge_cleaned'] = ".cleaned" if config['run_modes'][run_mode]['clean_dge'] else ""
    sample_config['run_mode']['polyA_adapter_trimmed'] = ".polyA_adapter_trimmed" if config['run_modes'][run_mode]['polyA_adapter_trimming'] else ""
    sample_config['run_mode']['mm_included'] = ".mm_included" if config['run_modes'][run_mode]['count_mm_reads'] else ""
    
    return sample_config

def get_all_tiles_dge_path(sample_config, sample_metadata, check_exists=False):
    puck_ids = sample_metadata['puck_barcode_file_id']
    puck_dges = []
    for puck_id in puck_ids:
        puck_dges.append(SMK_DGE_FILE.format(project_id=sample_metadata['project_id'],
                                             sample_id=sample_metadata['sample_id'],
                                             tile=f"_{puck_id}",
                                             **sample_config['run_mode']))
        check_file_exists(puck_dges[-1], exception=check_exists)

    return puck_dges, puck_ids

def get_stitched_dge(sample_config, sample_metadata, check_exists=False):
    _file = SMK_MULTIMODAL_PROCESSED_DGE.format(project_id=sample_metadata['project_id'],
                                        sample_id=sample_metadata['sample_id'])
    
    check_file_exists(_file, exception=check_exists)
    return _file

# TODO: this is a directory, check directory
def get_raw_image_dir_path(sample_config, sample_metadata, check_exists=False):
    _file = SMK_IMAGES_RAW.format(project_id=sample_metadata['project_id'],
                                  sample_id=sample_metadata['sample_id'])
    
    check_directory_exists(_file, exception=check_exists)
    return _file

def get_stitched_image_path(sample_config, sample_metadata, check_exists=False):
    _file = SMK_IMAGES_PROCESSED_STITCHED.format(project_id=sample_metadata['project_id'],
                                                 sample_id=sample_metadata['sample_id'])
    
    check_file_exists(_file, exception=check_exists)
    return _file

def _run_spatial_stitch(sample_config, sample_metadata):
    from openst.cli import get_spatial_stitch_parser
    puck_dges, puck_ids = get_all_tiles_dge_path(sample_config, sample_metadata, check_exists=True)
    required_arguments = ["--tiles"] + puck_dges + ["--tile-id"] + puck_ids

    return required_arguments


SUBPARSER_ACTIONS = {
    'image_stitch': None,
    'spatial_stitch': _run_spatial_stitch,
    'pairwise_aligner': None,
    'apply_transform': None,
    'manual_pairwise_aligner': None,
    'segment': None,
    'segment_merge': None,
    'transcript_assign': None,
    'merge_modalities': None,
    'pseudoimage': None,
    'preview': None
}

def _run_from_spacemake(parser, args, unknown_args):
    from openst.cli import setup_from_spacemake_subparsers
    
    if not os.path.isfile(config_path) or not os.path.isfile(project_df):
        raise SpacemakeError("""There is no 'config.yaml' and/or 'project_df.csv' in this directory.
                             Is it a valid spacemake directory? Run 'spacemake init --help' for more details.""")
    
    pdf = load_project_df()
    config = load_config()
    subcommand = unknown_args[0]

    sample_metadata = get_sample_metadata(pdf, args.project_id, args.sample_id)
    sample_config = get_config_from_sample(sample_metadata, config, args.run_mode)

    from_spacemake_subparsers = parser.add_subparsers(help="sub-command help", dest="subcommand")
    setup_from_spacemake_subparsers(from_spacemake_subparsers)

    # We get the arguments for from_spacemake
    arguments_list = []
    for key, value in args.__dict__.items():
        if key in ['subcommand', 'func']:
            continue
        arguments_list.extend([f"--{key.replace('_', '-')}", str(value)])
    
    # In the unknown_args, there will be the subcommand (because it was not added yet)
    required_arguments = SUBPARSER_ACTIONS[subcommand](sample_config, sample_metadata)
    args = parser.parse_args(arguments_list + unknown_args + required_arguments)

    logging.info(f"openst {args.subcommand} - running with the following parameters:")
    logging.info(args.__dict__)
    args.func(args)

if __name__ == "__main__":
    from openst.cli import get_from_spacemake_parser
    args = get_from_spacemake_parser().parse_args()
    _run_from_spacemake(args)
