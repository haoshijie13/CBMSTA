#!/usr/bin/env python3
from skimage import io
import numpy as np
import pandas as pd
from scipy.spatial import KDTree
import os
import argparse
import matplotlib.pyplot as plt
from scipy.spatial import distance
import scanpy as sc
import scipy.ndimage as nd

'''min = np.min(img)
max = np.max(img)
image_8bit = np.array(np.rint(255 * ((img - min) / (max - min))), dtype=np.uint8)'''


def main(args):
    if args.step != 'image' and args.step != 'distance':
        raise(ValueError("Please select either 'image' or 'distance' step."))
    elif args.step == 'image':
        #load data
        gem_data  = sc.read_h5ad(args.gem)
        gem_data = gem_data.obs.copy()
        gem_data = gem_data[['coor_x', 'coor_y', 'nCount_RNA']]
        gem_data = gem_data.astype('int')
        coords = np.zeros((gem_data['coor_y'].max() + 10, gem_data['coor_x'].max() + 10), dtype=int)
        exp_max = gem_data['nCount_RNA'].mean() * 2
        coords[gem_data['coor_y'], gem_data['coor_x']] = gem_data['nCount_RNA']
        coords[coords >= exp_max] = exp_max
        coords = coords * (255 / exp_max)
        coords = coords.astype('uint8')

        if not os.path.exists("mask_pre"):
            os.mkdir("mask_pre")
        io.imsave(f'mask_pre/{args.output}.bin50.tif', coords)



    elif args.step == 'distance':
        line = pd.read_csv(args.line,sep='\t', header=0, compression='infer', comment='#')
        line.columns = ['x', 'y', 'line_id']
        area_colour = np.loadtxt(args.colour).astype('int')



        gem_data = sc.read_h5ad(args.gem)
        gem_data = gem_data.obs.copy()
        gem_data = gem_data[['coor_x', 'coor_y', 'nCount_RNA']]
        gem_data = gem_data.astype('int')
        mask = io.imread(args.mask)
        # gem_data = pd.read_csv(args.gem,sep='\t', header=0, compression='infer', comment='#')
        # area_colour = args.colour

        if args.type == 'mouse':
            affineR = np.matrix(np.array([[0.5, 0, -0.5], [0, 0.5, -0.5], [0, 0, 1]]))
            mask = nd.affine_transform(mask.T, affineR, output_shape=[i * 2 for i in mask.T.shape], order=0)
            mask = mask.T
            mask = mask.astype('uint8')
            line['x'] = line['x']*2
            line['y'] = line['y']*2
        gem_cell = gem_data.copy()
        gem_cell['area'] = mask[gem_data['coor_y'], gem_data['coor_x']]
        gem_cell = gem_cell[gem_cell['area'] != 0]

        if not os.path.exists("distance_result"):
            os.mkdir("distance_result")


        all_area = pd.DataFrame(columns=['coor_x','coor_y','area','pdistance','plength','label'])
        for i in range(0, int(len(area_colour)/3)):
            # choice area coor and line
            m_label = area_colour[i*3]
            g_label = area_colour[i*3+1]
            w_label = area_colour[i*3+2]
            label_line = line[line['line_id'] == i+1].copy()
            m_area = gem_cell[gem_cell['area'] == m_label].copy()
            g_area = gem_cell[gem_cell['area'] == g_label].copy()
            w_area = gem_cell[gem_cell['area'] == w_label].copy()
            m_area['area'] = 'molecular'
            g_area['area'] = 'granular'
            w_area['area'] = 'white'

            label_line = np.array(label_line[['x', 'y']].copy(),dtype='float')
            m_coor = np.array(m_area[['coor_x', 'coor_y']].copy(),dtype='int')
            g_coor = np.array(g_area[['coor_x', 'coor_y']].copy(), dtype='int')
            w_coor = np.array(w_area[['coor_x', 'coor_y']].copy(), dtype='int')

            length = np.zeros(len(label_line))
            for j in range(1, len(label_line)):
                length[j] = distance.euclidean(label_line[j, :], label_line[j - 1, :])
            length = np.cumsum(length)

            # find distance
            kdtree = KDTree(label_line)
            darraym, iarraym = kdtree.query(m_coor)
            darrayg, iarrayg = kdtree.query(g_coor)
            darrayw, iarrayw = kdtree.query(w_coor)

            m_area['pdistance'] = darraym
            g_area['pdistance'] = darrayg
            w_area['pdistance'] = darrayw

            m_area['plength'] = length[np.array(iarraym)]
            g_area['plength'] = length[np.array(iarrayg)]
            w_area['plength'] = length[np.array(iarrayw)]

            m_area['label'] = i+1
            g_area['label'] = i+1
            w_area['label'] = i+1

            all_area = pd.concat([all_area,m_area,g_area,w_area])

        all_area.to_csv(f'distance_result/{args.output}.raw.distance.csv', sep='\t', header=True, index=False)
        img1 = all_area.plot.scatter(x='coor_x', y='coor_y', c='pdistance', colormap='viridis')
        plt.savefig(f'distance_result/{args.output}.raw1.result.png')
        plt.close()

        distance1 = all_area.copy()
        distance1.loc[(distance1['area'] == 'granular') | (distance1['area'] == 'white'), 'pdistance'] = distance1[(distance1['area'] == 'granular') | (distance1['area'] == 'white')]['pdistance'].copy() * -1
        img0 = distance1.plot.scatter(x='coor_x', y='coor_y', c='pdistance', colormap='viridis')
        plt.savefig(f'distance_result/{args.output}.raw2.result.png')
        plt.close()

        plot = all_area[['coor_x','coor_y','area','pdistance','plength','label']].copy()
        #plot['pdistance'] = (plot['pdistance']+0.5).astype(int)
        plot['plength'] = (plot['plength']+0.5).astype(int)
        new = plot.groupby(['plength', 'label','area']).agg(distancemin=('pdistance', 'min'),distancemax=('pdistance', 'max')).reset_index()
        n_plot = pd.merge(plot,new)


        #p_max = plot['pdistance'].max()
        molecular = n_plot[n_plot['area'] == 'molecular'].copy()
        white = n_plot[n_plot['area'] == 'white'].copy()
        granular = n_plot[n_plot['area'] == 'granular'].copy()
        #m_max = molecular['pdistance'].max()
        w_max = white['pdistance'].max()
        w_min = white['pdistance'].min()
        #g_max = granular['pdistance'].max()

        molecular['pdistance'] = (molecular['pdistance']-molecular['distancemin'])/(molecular['distancemax']-molecular['distancemin']+0.000001)
        #white['pdistance'] = (((white['pdistance']-white['distancemin'])/(white['distancemax']-white['distancemin']+0.000001))+1) *-1
        white['pdistance'] = ((white['pdistance']-w_min)/(w_max-w_min)+1)*-1
        granular['pdistance'] = ((granular['pdistance']-granular['distancemin'])/(granular['distancemax']-granular['distancemin']+0.000001)) *-1
        norm = pd.concat([molecular,white,granular])
        norm = norm.fillna(0)
        norm.to_csv(f'distance_result/{args.output}.norm.distance.csv', sep='\t', header=True, index=False)
        img = norm.plot.scatter(x='coor_x', y='coor_y', c='pdistance', colormap='viridis')
        plt.savefig(f'distance_result/{args.output}.norm.result.png')
        plt.close()
        print('mission complete')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--step", help="either 'image' or 'distance' ", type=str, default="distance")
    parser.add_argument("-t", "--type", help="either 'other' or 'mouse' ", type=str, default="other")
    parser.add_argument("-m", "--mask", help="input the mask image", type=str, default='', nargs='?')
    parser.add_argument("-g", "--gem", help="input the gem data", type=str, default='', nargs='?')
    parser.add_argument("-c", "--colour", help="input the colour list, follow molecular1, granular1, white1, molecular2...,if non-exist input 0", type=str, default='', nargs='?')
    parser.add_argument("-l", "--line", help="input pos file path", type=str, default='', nargs='?')
    parser.add_argument("-o", "--output", help="directory to save files", type=str, default='result')
    args = parser.parse_args()
    main(args)

