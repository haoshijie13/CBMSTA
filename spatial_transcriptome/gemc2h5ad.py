#!/usr/bin/env python3
import argparse
import os

def main(args):
    from scipy.sparse import csr_matrix
    from anndata import AnnData
    import pandas as pd

    gem_data = pd.read_csv(args.gem, sep='\s+', header=0, compression='infer', comment='#',low_memory=True,na_filter=False,engine='c')
    gem_data.columns = args.columns
    gem_data = gem_data[gem_data['cell']!=0].copy()
    coor_cell = gem_data.groupby(by=['cell'])['x','y'].mean().reset_index()
    coor_cell['cell_coor'] = coor_cell['x'].map(str) + '_' + coor_cell['y'].map(str)

    gem_data = gem_data[['geneID', 'x', 'y', 'MIDCount', 'cell']]
    gem_data['cell'] = gem_data['cell'].map(dict(zip(coor_cell['cell'], coor_cell['cell_coor'])))
    gem_data = gem_data.groupby(['cell', 'geneID']).agg(MIDCount=('MIDCount', 'sum')).reset_index()


    uniq_cell, uniq_gene = gem_data.cell.unique(), gem_data.geneID.unique()
    uniq_cell, uniq_gene = list(uniq_cell), list(uniq_gene)
    cell_dict = dict(zip(uniq_cell, range(0, len(uniq_cell))))
    gene_dict = dict(zip(uniq_gene, range(0, len(uniq_gene))))
    gem_data["csr_x_ind"] = gem_data["cell"].map(cell_dict)
    gem_data["csr_y_ind"] = gem_data["geneID"].map(gene_dict)

    matrix = csr_matrix((gem_data['MIDCount'], (gem_data["csr_x_ind"], gem_data["csr_y_ind"])),
                        shape=((len(uniq_cell), len(uniq_gene))))

    var = pd.DataFrame({"gene_short_name": uniq_gene})
    var.set_index("gene_short_name", inplace=True)
    obs = pd.DataFrame({"cell": gem_data['cell'].unique().tolist()})
    obs["coor_x"] = pd.Series(obs['cell']).str.split('_', expand=True).iloc[:, 0]
    obs['coor_y'] = pd.Series(obs['cell']).str.split('_', expand=True).iloc[:, 1]
    obs['cell'] = obs['cell'].map(dict(zip(coor_cell['cell_coor'], coor_cell['cell'])))
    obsm = {"spatial": obs.loc[:, ['coor_x', "coor_y"]].values}
    obs.set_index("cell", inplace=True)

    adata = AnnData(matrix, dtype=matrix.dtype, obs=obs.copy(), var=var.copy(), obsm=obsm.copy())

    if not os.path.exists("STh5ad"):
        print("create STh5ad")
        os.mkdir("STh5ad")
    adata.write_h5ad(f'STh5ad/{args.output}.h5ad', compression='lzf')
    print(f'file save in STh5ad/{args.output}.h5ad')
    print('mission completed')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gem", help="input the gem file path", type=str, default='')
    parser.add_argument("-o", "--output", help="output name, default output", type=str, default='output',nargs='?')
    parser.add_argument("-C", "--columns",
                        help="set the gem columns names, must include(geneID x y MIDCount cell) default: geneID x y MIDCount cell",
                        type=str, default=['geneID', 'x', 'y', 'MIDCount', 'cell'], nargs='+')
    args = parser.parse_args()
    main(args)
