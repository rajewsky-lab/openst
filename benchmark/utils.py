import os
import pysam
from tqdm import tqdm
import numpy as np
import subprocess

def downsampling(base_path, max_reads, steps, override=False):
    os.makedirs(os.path.dirname(base_path + 'bam_downsampling/'), exist_ok=True)
    subsets_n_reads = np.linspace(max_reads/steps, max_reads, steps)
    reads = []
    with pysam.AlignmentFile(base_path + f'fov_subset_1mm2_area_reads_collated.bam', 'rb') as f:
        for i, read in tqdm(enumerate(f)):
            if i > max_reads:
                break
            reads.append(read)
            if i in subsets_n_reads:
                print(f'processed {i} reads')
                if os.path.exists(base_path + f'bam_downsampling/fov_subset_1mm2_area_reads_collated_{i // 1000000}M.bam') and not override:
                    print('file already exists')
                    continue
                subset_file = pysam.AlignmentFile(base_path + f'bam_downsampling/fov_subset_1mm2_area_reads_collated_{i // 1000000}M.bam', 'wb', template=f)
                for r in reads:
                    subset_file.write(r)
                subset_file.close()
            
    return
    

def digital_expression(base_path, subset, override=False, with_intronic=True):
    os.makedirs(os.path.dirname(base_path + 'downsampling_results/'), exist_ok=True)
    if with_intronic:
        os.makedirs(os.path.dirname(base_path + 'downsampling_results/with_intronic/'), exist_ok=True)
        
    n_subset = subset.split('_')[-1].split('M')[0]
    input_path = base_path + 'bam_downsampling/' + subset
    cell_bc_path = base_path + 'fov_subset_1mm2_area_bcs.txt'
    
    if with_intronic:
        output_path = base_path + 'downsampling_results/with_intronic/' + f'dropseqtools_dge_{n_subset}M.txt.gz'
        summary_path = base_path + 'downsampling_results/with_intronic/' + f'dropseqtools_summary_{n_subset}M.txt'
        if os.path.exists(output_path) and not override:
            print(f'subset {n_subset}M has already been processed')
            return
        print(f'processing {n_subset}M subset')
        
        cmd = ['./../../Drop-seq_tools-2.5.1/DigitalExpression',
                f'--INPUT {input_path}',
                f'--OUTPUT {output_path}',
                f'--CELL_BC_FILE {cell_bc_path}',
                f'--SUMMARY {summary_path}',
                '--CELL_BARCODE_TAG CB',
                '--MOLECULAR_BARCODE_TAG MI',
                '--TMP_DIR /tmp',
                '--LOCUS_FUNCTION_LIST CODING',
                '--LOCUS_FUNCTION_LIST UTR',
                '--LOCUS_FUNCTION_LIST INTRONIC'
                ]
        
    else:
        output_path = base_path + 'downsampling_results/' + f'dropseqtools_dge_{n_subset}M.txt.gz'
        summary_path = base_path + 'downsampling_results/' + f'dropseqtools_summary_{n_subset}M.txt'
        if os.path.exists(output_path) and not override:
            print(f'subset {n_subset}M has already been processed')
            return
        print(f'processing {n_subset}M subset')
        
        cmd = ['./../../Drop-seq_tools-2.5.1/DigitalExpression',
                f'--INPUT {input_path}',
                f'--OUTPUT {output_path}',
                f'--CELL_BC_FILE {cell_bc_path}',
                f'--SUMMARY {summary_path}',
                '--CELL_BARCODE_TAG CB',
                '--MOLECULAR_BARCODE_TAG MI',
                '--TMP_DIR /tmp']   

    subprocess.call(cmd)
        
        
        
def get_genic_reads_and_umis(base_path, max_reads, steps):
    subsets_n_reads = np.linspace(max_reads/steps, max_reads, steps)
    total_counts = []
    total_genic_reads = []
    files = os.listdir(base_path)
    n_files = len(files) // 2
    for file in files:
        n_subset = file.split('_')[-1].split('M')[0]
        if int(n_subset) * 1000000 not in subsets_n_reads or 'summary' not in file:
            continue
        
        with open(base_path + f'/dropseqtools_summary_{n_subset}M.txt', 'r') as f:
            lines = f.readlines()[7:]
            genic_reads_subset = 0
            counts_subset = 0
            for l in lines:
                l_split = l.split('\t')
                if len(l_split) != 4:
                    print(l)
                    continue
                bc, genic_reads, counts, genes = l_split
                genic_reads_subset += int(genic_reads)
                counts_subset += int(counts)

            total_counts.append(counts_subset)
            total_genic_reads.append(genic_reads_subset)
    
    total_counts = sorted(total_counts)
    total_genic_reads = sorted(total_genic_reads)
        
    return total_counts, total_genic_reads