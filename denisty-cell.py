#!/usr/bin/env python3
import numpy as np
import pandas as pd
import sys
from scipy import stats
#from collections import Counter

distance = sys.argv[1]
cell_data = sys.argv[2]
slice_name = sys.argv[3]

distance = pd.read_csv(distance,sep="\t",header=0)
distance = distance[(distance['area']=='granular')|(distance['area']=='molecular')]
distance['coor'] = distance['coor_x'].map(str)+'_'+ distance['coor_y'].map(str)

cell_data = pd.read_csv(cell_data,sep=",",header=0)
cell_data['x'] = round(cell_data['x']/50,0).astype(int)
cell_data['y'] = round(cell_data['y']/50,0).astype(int)
cell_data['coor'] = cell_data['x'].map(str)+'_'+ cell_data['y'].map(str)

cell_data['distance'] = cell_data['coor'].map(dict(zip(distance['coor'],distance['pdistance'])))
cell_data['cell'] = cell_data['celltype_pred_filter']

cell = ["Fibroblast", "MLI1", "MLI2", "Bergmann", "Purkinje","PLI" ,"Granule", "Golgi", "UBC", "OPC", "Astrocyte", "ODC"]
#col = ['distance']+cell
#col_index = [i in col for i in obs]
#if not all(col_index):
#    print(f"{slice_name} has't distance or cell")
#    quit()


print("adata has been loaded")
all = cell_data[['distance','cell']].copy()
all = all.dropna()
all = all.reset_index(drop=True)
#all = pd.merge(distance, bin, on=['coor_x', 'coor_y'])

X_impuse = np.linspace(-1,1,201)
slice = np.repeat(slice_name, len(X_impuse))
slice_df = pd.DataFrame({'slice':slice})
all_matrix = [X_impuse]
#all['num'] = 1
#all['distance'] = all['distance']*100
#all['distance'] = all['distance'].astype(int)
#all = all.groupby(["distance"])[cell+['num']].sum().reset_index()
#all['distance'] = all['distance']/100
#all = all.fillna(0)

print(f'{slice_name} start')
for i in range(len(cell)):
    name = cell[i]
    print(name)
    #expression_gene = (all[name]+0.001).tolist()
    tmp = all[all['cell']==name].copy()
    kde = stats.gaussian_kde(np.array(tmp["distance"]))
    density = kde(X_impuse)

    #normalize
    bw = kde.scotts_factor()
    num_kde = stats.gaussian_kde(np.array(all["distance"]),bw_method=bw)
    #num_kde = stats.gaussian_kde(np.array(all["distance"]), bw_method=bw)
    #print(f'{name}:{bw}')
    num_density = num_kde(X_impuse)
    all_matrix = all_matrix + [density/num_density]

n_slice_df = pd.DataFrame(np.array(all_matrix).T,columns=['distance']+cell)
n_slice_df = pd.concat([slice_df,n_slice_df],axis=1)
print(f"{slice_name} complete")

n_slice_df.to_csv(f"{slice_name}_cell_density.csv", sep='\t', header=True, index=False)
print("mission complete")

