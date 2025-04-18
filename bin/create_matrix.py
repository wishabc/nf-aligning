import pandas as pd
import sys
import numpy as np
import scipy.sparse as sp
from tqdm import tqdm


def main(samples_order, filelist_map):
    print("Reading files...")
    f = np.load(filelist_map[samples_order[0]])
    data = np.zeros((len(samples_order), f.shape[0]), dtype=np.float16)
    for i, sample_id in enumerate(tqdm(samples_order)):
        file = filelist_map[sample_id]
        try:
            f = np.load(file).astype(np.float16)
            data[i, :] = f
        except:
            print("Problems with", sample_id)
            raise
    
    total_els = data.shape[0] * data.shape[1]
    print(total_els, (data > 0).sum() / total_els)
    # data = sp.coo_matrix(data).tocsr()
    # sparsity = 1.0 - data.nnz / (total_els)
    # print(f"Sparsity: {sparsity:.2%}")
    return data


if __name__ == '__main__':
    samples_order = pd.read_table(sys.argv[1])['ag_id'].values
    filelist = np.loadtxt(sys.argv[2], dtype=str)
    ids_map = {x.replace('.max_pvals.npy', ''): x for x in filelist}
    sparse_res = main(samples_order, ids_map)
    np.save(sys.argv[3], sparse_res)