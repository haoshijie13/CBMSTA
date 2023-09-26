#!/usr/bin/env python3
import numpy as np
import pandas as pd
import sys
from scipy import stats
#from collections import Counter

species = sys.argv[1]
adata = sys.argv[2]
slice_name = sys.argv[3]

print("script start")
#distance = pd.read_csv(distance, sep='\t', header=0, compression='infer', comment='#')
#gene_raw_df = pd.read_csv("/jdfssz3/ST_STOMICS/P22Z10200N0661/huangzhi/cerebellum/purkinje_d/paper/mart_export.humanMacaqeMarmosetMouse.oneToOneOrth.ensembl91.20220428.txt", sep="\t", header=0)
#gene_raw_df = pd.read_csv("/jdfssz3/ST_STOMICS/P22Z10200N0661/huangzhi/cerebellum/purkinje_d/paper/test.txt", sep="\t", header=0)
#gene_raw_list = gene_raw_df[f'{species}Gene'].tolist()
# cluster_name_list = args.cluster.copy()
import scanpy as sc
adata = sc.read_h5ad(adata)
if 'distance' not in adata.obs.columns:
    print(f"{slice_name} has't distance")
    quit()
adata = adata[(adata.obs['area']=='granular')|(adata.obs['area']=='molecular')]
obs = adata.obs.copy()
obs.index = list(range(0, len(obs)))
adata.var = adata.var.set_index('features')
adata.var.index = adata.var.index.astype('str')
print("adata has been loaded")

#transform_geneid = adata.var.index.map(dict(zip(gene_raw_df[f'{species}Gene'],gene_raw_df['marmosetGene'])))
gene_expression = np.array(adata.X.sum(axis=0))[0]
gene_index = gene_expression > 100
gene_matrix = adata.X[:, gene_index]
gene_list = adata.var[gene_index].index.tolist()
df_matrix = pd.DataFrame(gene_matrix.toarray(), columns=gene_list)
obs = obs[["distance","area"]].copy()
all = pd.concat([obs, df_matrix], axis=1)
#all = pd.merge(distance, bin, on=['coor_x', 'coor_y'])

save = pd.DataFrame({"gene":adata.var.index,"expression":gene_expression})
save['slice'] = slice_name
save.to_csv(f'{slice_name}.merge.csv', sep='\t', header=True, index=True)
print("merge file has been saved")

X_impuse = np.linspace(-1,1,201)
slice = np.repeat(slice_name, len(X_impuse))
slice_df = pd.DataFrame({'slice':slice})
all_matrix = [X_impuse]
all['num'] = 1
all['distance'] = all['distance']*100
all['distance'] = all['distance'].astype(int)
all = all.groupby(["distance"])[gene_list+['num']].sum().reset_index()
all['distance'] = all['distance']/100
#all = all.fillna(0)

print(f'{slice_name} start')
for i in range(len(gene_list)):
    name = gene_list[i]
    #if len(all[all[name]!=0])==1:
        #print(f'{name} only in one distance')
        #all_matrix = all_matrix + [np.repeat(0,len(X_impuse))]
        #continue
    expression_gene = (all[name]+0.001).tolist()
    kde = stats.gaussian_kde(np.array(all["distance"]),weights=expression_gene)
    density = kde(X_impuse)

    #normalize
    bw = kde.scotts_factor()
    num_kde = stats.gaussian_kde(np.array(all["distance"]), weights=all['num'].tolist(),bw_method=bw)
    #num_kde = stats.gaussian_kde(np.array(all["distance"]), bw_method=bw)
    #print(f'{name}:{bw}')
    num_density = num_kde(X_impuse)
    all_matrix = all_matrix + [density/num_density]

n_slice_df = pd.DataFrame(np.array(all_matrix).T,columns=['distance']+gene_list)
n_slice_df = pd.concat([slice_df,n_slice_df],axis=1)
print(f"{slice_name} complete")

n_slice_df.to_csv(f"{slice_name}_density.csv", sep='\t', header=True, index=False)
print("mission complete")

