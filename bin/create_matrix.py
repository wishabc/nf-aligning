import pandas as pd
import sys
import numpy as np
import scipy.sparse as sp
from tqdm import tqdm

def main(samples_order, filelist_map):
    data = []
    print("Reading files...")
    for sample_id in tqdm(samples_order):
        file = filelist_map[sample_id]
        df = pd.read_table(file, names=['chrom', 'start', 'end', 'neglog10_pval'], header=None)
        data.append(df['neglog10_pval'].values.astype(np.float32))
    
    data = np.vstack(data)
    data = sp.coo_matrix(data).tocsr()
    sparsity = 1.0 - data.nnz / (data.shape[0] * data.shape[1])
    print(f"Sparsity: {sparsity:.2%}")
    return data


if __name__ == '__main__':
    samples_order = pd.read_table(sys.argv[1])['ag_id'].values
    filelist = np.loadtxt(sys.argv[2], dtype=str)
    ids_map = {x.replace('.max_pvals.bed', ''): x for x in filelist}
    sparse_res = main(samples_order, ids_map)
    sp.save_npz(sys.argv[3], sparse_res)