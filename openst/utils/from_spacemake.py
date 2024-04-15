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
SMK_RAW = 'raw_data'

SMK_DGE = os.path.join(SMK_PROJECT, SMK_SAMPLE_PROCESSED, 'illumina/complete_data/dge')
SMK_DGE_FILE = os.path.join(SMK_DGE, 'dge{dge_type}{dge_cleaned}{polyA_adapter_trimmed}{mm_included}.spatial_beads{tile}.h5ad')
SMK_STITCHED_DGE = os.path.join(SMK_DGE, 'dge{dge_type}{dge_cleaned}{polyA_adapter_trimmed}{mm_included}.spatial_beads_puck_collection.h5ad')

SMK_IMAGES_PROCESSED = os.path.join(SMK_PROJECT, SMK_SAMPLE_PROCESSED, 'images')
SMK_IMAGES_PROCESSED_STITCHED = os.path.join(SMK_IMAGES_PROCESSED, 'Image_Stitched_Composite.tif')
SMK_IMAGES_RAW = os.path.join(SMK_PROJECT, SMK_RAW, 'images')

SMK_MULTIMODAL_PROCESSED = os.path.join(SMK_PROJECT, SMK_SAMPLE_PROCESSED, 'multimodal')
SMK_MULTIMODAL_PROCESSED_DGE = os.path.join(SMK_MULTIMODAL_PROCESSED, 'stitched_spots.h5ad')
SMK_MULTIMODAL_PROCESSED_DGE_SEGMENTED = os.path.join(SMK_MULTIMODAL_PROCESSED, 'stitched_segmented.h5ad')

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
    elif len(_sample) == 0:
        raise SpacemakeError(f"Entry project_id='{project_id}' and sample_id='{sample_id}' not found")
    
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

def get_key_or_default(config, category='run_modes', name='default', key='count_intronic_reads'):
    if key in config[category][name].keys():
        return config[category][name][key]
    else:
        return config[category]['default'][key]

def get_config_from_sample(sample_metadata, config, run_mode=""):
    sample_config = {}

    # Puck
    sample_config['puck'] = config['pucks'][sample_metadata['puck']]

    # Run mode options
    if run_mode is "" and len(sample_metadata['run_mode']) > 1:
        raise SpacemakeError(f"""There are various 'run_modes' for this sample:
                             \t{sample_metadata['run_mode']}
                              Specify only one using the '--run-mode' argument""")
    elif run_mode is "":
        run_mode = sample_metadata['run_mode'][0]

    sample_config['run_mode'] = {}
    sample_config['run_mode']['dge_type'] = ".all" if get_key_or_default(config, 'run_modes', run_mode, 'count_intronic_reads') else ".exon"
    sample_config['run_mode']['dge_cleaned'] = ".cleaned" if get_key_or_default(config, 'run_modes', run_mode, 'clean_dge') else ""
    sample_config['run_mode']['polyA_adapter_trimmed'] = ".polyA_adapter_trimmed" if get_key_or_default(config, 'run_modes', run_mode, 'polyA_adapter_trimming') else ""
    sample_config['run_mode']['mm_included'] = ".mm_included" if get_key_or_default(config, 'run_modes', run_mode, 'count_mm_reads') else ""
    
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


def get_multimodal_dge(sample_config, sample_metadata, check_exists=False):
    _file = SMK_MULTIMODAL_PROCESSED_DGE.format(project_id=sample_metadata['project_id'],
                                        sample_id=sample_metadata['sample_id'])
    
    if not check_directory_exists(_file):
        _parent_dir = os.path.dirname(_file)
        logging.info(f"Creating directory at {_parent_dir}")
        os.mkdir(_parent_dir)

    check_file_exists(_file, exception=check_exists)
    return _file

def get_stitched_dge(sample_config, sample_metadata, check_exists=False):
    _file = SMK_STITCHED_DGE.format(project_id=sample_metadata['project_id'],
                                    sample_id=sample_metadata['sample_id'],
                                    **sample_config['run_mode'])

    if not check_directory_exists(_file):
        _parent_dir = os.path.dirname(_file)
        logging.info(f"Creating directory at {_parent_dir}")
        os.mkdir(_parent_dir)

    check_file_exists(_file, exception=check_exists)
    return _file

def get_stitched_segmented_dge(sample_config, sample_metadata, check_exists=False):
    _file = SMK_MULTIMODAL_PROCESSED_DGE_SEGMENTED.format(project_id=sample_metadata['project_id'],
                                        sample_id=sample_metadata['sample_id'])
    
    if not check_directory_exists(_file):
        _parent_dir = os.path.dirname(_file)
        logging.info(f"Creating directory at {_parent_dir}")
        os.mkdir(_parent_dir)

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
    
    if not check_directory_exists(_file):
        _parent_dir = os.path.dirname(_file)
        logging.info(f"Creating directory at {_parent_dir}")
        os.mkdir(_parent_dir)

    check_file_exists(_file, exception=check_exists)
    return _file

def _run_spatial_stitch(sample_config, sample_metadata):
    SMK_STITCHED_DGE

    puck_dges, puck_ids = get_all_tiles_dge_path(sample_config, sample_metadata, check_exists=True)
    required_arguments = ["--tiles"] + puck_dges + ["--tile-id"] + puck_ids
    output_file = get_stitched_dge(sample_config, sample_metadata)
    required_arguments += ["--h5-out", output_file]

    return required_arguments

def _run_image_stitch(sample_config, sample_metadata):
    raw_image_dir = get_raw_image_dir_path(sample_config, sample_metadata, check_exists=True)
    image_out = get_stitched_image_path(sample_config, sample_metadata)
    required_arguments = ["--image-indir", raw_image_dir, "--image-out", image_out]

    return required_arguments

def _run_segment_merge(sample_config, sample_metadata):
    stitched_dge = get_multimodal_dge(sample_config, sample_metadata, check_exists=True)
    required_arguments = ["--h5-in", stitched_dge]
    required_arguments += ["--mask-out", 'uns/spatial/staining_image_mask_merged']

    return required_arguments

def _run_segment(sample_config, sample_metadata):
    stitched_dge = get_multimodal_dge(sample_config, sample_metadata, check_exists=True)
    required_arguments = ["--image-in", 'uns/spatial/staining_image', "--h5-in"] + [stitched_dge]
    required_arguments += ["--mask-out", 'uns/spatial/staining_image_mask']

    return required_arguments

def _run_transcript_assign(sample_config, sample_metadata):
    stitched_dge = get_multimodal_dge(sample_config, sample_metadata, check_exists=True)
    stitched_segmented_dge = get_stitched_segmented_dge(sample_config, sample_metadata)
    required_arguments = ["--mask-in", 'uns/spatial/staining_image_mask', "--h5-in"] + [stitched_dge]
    required_arguments += ["--h5-out", stitched_segmented_dge]

    return required_arguments

def _run_merge_modalities(sample_config, sample_metadata):
    import shutil

    stitched_dge = get_stitched_dge(sample_config, sample_metadata, check_exists=True)
    image_out = get_stitched_image_path(sample_config, sample_metadata, check_exists=True)
    merged_dge = get_multimodal_dge(sample_config, sample_metadata)

    if not check_directory_exists(merged_dge):
        _parent_dir = os.path.dirname(merged_dge)
        logging.info(f"Creating directory at {_parent_dir}")
        os.mkdir(_parent_dir)

    if check_file_exists(merged_dge, exception=False):
        logging.warn(f"No need to create {merged_dge} - it exists")  
    else:
        logging.info(f"Copying {stitched_dge} into {merged_dge}")  
        shutil.copy(stitched_dge, merged_dge)
    
    required_arguments = ["--h5-in", merged_dge, "--image-in", image_out]

    return required_arguments

def _run_manual_pairwise_aligner(sample_config, sample_metadata):
    stitched_dge = get_multimodal_dge(sample_config, sample_metadata, check_exists=True)
    required_arguments = ["--h5-in", stitched_dge]

    return required_arguments

def _run_apply_transform(sample_config, sample_metadata):
    stitched_dge = get_multimodal_dge(sample_config, sample_metadata, check_exists=True)
    required_arguments = ["--h5-in", stitched_dge]

    return required_arguments

def _run_pseudoimage(sample_config, sample_metadata):
    stitched_dge = get_multimodal_dge(sample_config, sample_metadata, check_exists=True)
    required_arguments = ["--h5-in", stitched_dge]

    return required_arguments

def _run_preview(sample_config, sample_metadata):
    stitched_dge = get_multimodal_dge(sample_config, sample_metadata, check_exists=True)
    required_arguments = ["--h5-in", stitched_dge]

    return required_arguments

def _run_pairwise_aligner(sample_config, sample_metadata):
    stitched_dge = get_multimodal_dge(sample_config, sample_metadata, check_exists=True)
    required_arguments = ["--h5-in", stitched_dge]

    return required_arguments


SUBPARSER_ACTIONS = {
    'image_stitch': _run_image_stitch,
    'spatial_stitch': _run_spatial_stitch,
    'pairwise_aligner': _run_pairwise_aligner,
    'apply_transform': _run_apply_transform,
    'manual_pairwise_aligner': _run_manual_pairwise_aligner,
    'segment': _run_segment,
    'segment_merge': _run_segment_merge,
    'transcript_assign': _run_transcript_assign,
    'merge_modalities': _run_merge_modalities,
    'pseudoimage': _run_pseudoimage,
    'preview': _run_preview
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

    if subcommand not in SUBPARSER_ACTIONS.keys():
        logging.error(f"Could not find subcommand '{subcommand}'. Please choose one of: {SUBPARSER_ACTIONS.keys()}")
        exit(1)
    
    # In the unknown_args, there will be the subcommand (because it was not added yet)
    required_arguments = SUBPARSER_ACTIONS[subcommand](sample_config, sample_metadata)
    args = parser.parse_args(arguments_list + [subcommand] + required_arguments + unknown_args[1:])

    logging.info(f"openst {args.subcommand} - running with the following parameters:")
    logging.info(args.__dict__)
    args.func(args)

if __name__ == "__main__":
    from openst.cli import get_from_spacemake_parser
    args = get_from_spacemake_parser().parse_args()
    _run_from_spacemake(args)
