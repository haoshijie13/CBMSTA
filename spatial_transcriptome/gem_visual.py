#!/usr/bin/env python3
import os
import argparse


def main(args):
    import pandas as pd
    import numpy as np
    import cv2
    print('script start')
    global gem_data, save_name

    if len(args.normalize) % 2 != 0:
        print('Error: you should imput even normalize')
        quit()
    if args.function != 'split' and args.function != 'all' and args.function != 'split' and args.function != 'nCount':
        print('please set -f all/merge/split/nCount')

    if args.function == 'nCount':
        gem_data = pd.read_csv(f'{args.gem}', sep='\s+', header=0, compression='infer',comment='#', low_memory=True, na_filter=False, engine='c')
        print('load data finish')

    else:
        all_list = args.red + args.blue + args.green
        grep_list = '|'.join(all_list)
        grep_num = np.random.randint(0,10000)
        os.system(f'zgrep -w -E "{grep_list}" {args.gem} > {args.output}.{grep_num}.gem_mask.txt')
        gem_data = pd.read_csv(f'{args.output}.{grep_num}.gem_mask.txt', sep='\s+', header=0, compression='infer', comment='#',low_memory=True,na_filter=False,engine='c')
        os.system(f'rm {args.output}.{grep_num}.gem_mask.txt')
        print('load data finish')

    gem_data.columns = args.columns
    gem_data = gem_data[['variable', 'x', 'y', 'value']].copy()
    if len(args.chop) != 0:
        gem_data = gem_data[(gem_data['x'] > args.chop[0]) & (gem_data['x'] < args.chop[1]) & (gem_data['y'] > args.chop[2]) & (gem_data['y'] < args.chop[3])].copy()
        gem_data['x'] = gem_data['x'] - gem_data['x'].min()
        gem_data['y'] = gem_data['y'] - gem_data['y'].min()

    gem_data['x'] = round(gem_data['x']/args.bin, 0).astype(int)
    gem_data['y'] = round(gem_data['y']/args.bin, 0).astype(int)
    gem_data = gem_data.groupby(['variable','x', 'y']).agg(value=('value', 'sum')).reset_index()
    if not os.path.exists("STimg"):
        os.mkdir("STimg")
    print('plot start, save figure in STimg')
    save_name = args.output


    if args.function == 'nCount':

        print("max_value: ",gem_data['value'].max()," min_value: ",gem_data['value'].min())
        if len(args.normalize)!=0:
            for i in range(0,int(len(args.normalize)/2)):
                coords = np.zeros((gem_data['y'].max()+ 10, gem_data['x'].max() + 10), dtype=float)
                coords[gem_data['y'], gem_data['x']] = gem_data['value']
                coords[coords <= args.normalize[2*i]] = args.normalize[2*i]
                coords[coords >= args.normalize[2*i+1]] = args.normalize[2*i+1]
                #coords = (coords - args.normalize[2*i])  * (255 / (args.normalize[2*i+1]-args.normalize[2*i]))
                coords = (coords - np.min(coords))*(255/(np.max(coords)-np.min(coords)))
                coords = coords.astype('uint8')
                cv2.imwrite(f'STimg/{save_name}.{args.normalize[2*i]}-{args.normalize[2*i+1]}.nCount.tif', coords)
        else:
            coords = np.zeros((gem_data['y'].max() + 10, gem_data['x'].max() + 10), dtype=float)
            coords[gem_data['y'], gem_data['x']] = gem_data['value']
            coords = (coords - np.min(coords)) * (255 / (np.max(coords) - np.min(coords)))
            coords = coords.astype('uint8')
            cv2.imwrite(f'STimg/{save_name}.nCount.tif', coords)
        print('nCount plot finish')


    if args.function == 'merge' or args.function == 'all':
        coords = np.zeros((gem_data['y'].max() + 10, gem_data['x'].max() + 10, 3), dtype=float)
        if len(args.red)!=0:
            red = gem_data[gem_data['variable'].isin(args.red)].copy()
            red = red.groupby(['x', 'y']).agg(value=('value', 'sum')).reset_index()
            rm = red['value'].mean() * args.contrast[0]
            coords[red['y'], red['x'], 2] = red['value']
            coords[coords[:, :, 2] >= rm,2] = rm
            coords[:, :, 2] = coords[:, :, 2] * (255 / np.max(coords[:, :, 2]))

        if len(args.green)!=0:
            green = gem_data[gem_data['variable'].isin(args.green)].copy()
            green = green.groupby(['x', 'y']).agg(value=('value', 'sum')).reset_index()
            gm = green['value'].mean() * args.contrast[1]
            coords[green['y'], green['x'], 1] = green['value']
            coords[coords[:, :, 1] >= gm,1] = gm
            coords[:, :, 1] = coords[:, :, 1] * (255 / np.max(coords[:, :, 1]))

        if len(args.blue)!=0:
            blue = gem_data[gem_data['variable'].isin(args.blue)].copy()
            blue = blue.groupby(['x', 'y']).agg(value=('value', 'sum')).reset_index()
            bm = blue['value'].mean() * args.contrast[2]
            coords[blue['y'], blue['x'], 0] = blue['value']
            coords[coords[:, :, 0] >= bm,0] = bm
            coords[:, :, 0] = coords[:, :, 0] * (255 / np.max(coords[:, :, 0]))

        coords = coords.astype('uint8')

        cv2.imwrite(f'STimg/{save_name}.merge.tif', coords)
        print('merge plot finish')

    if args.function == 'split' or args.function == 'all':

        if len(args.red) != 0:
            for i in args.red:
                coords = np.zeros((gem_data['y'].max() + 10, gem_data['x'].max() + 10, 3), dtype=float)
                red = gem_data[gem_data['variable']==i].copy()
                rm = red['value'].mean() * args.contrast[0]
                coords[red['y'], red['x'], 2] = red['value']
                coords[coords[:, :, 2] >= rm,2] = rm
                coords[:, :, 2] = coords[:, :, 2] * (255 / np.max(coords[:, :, 2]))
                coords = coords.astype('uint8')
                cv2.imwrite(f'STimg/{save_name}.{i}.red.tif', coords)

        if len(args.green) != 0:
            for i in args.green:
                coords = np.zeros((gem_data['y'].max() + 10, gem_data['x'].max() + 10, 3), dtype=float)
                green = gem_data[gem_data['variable']==i].copy()
                gm = green['value'].mean() * args.contrast[1]
                coords[green['y'], green['x'], 1] = green['value']
                coords[coords[:, :, 1] >= gm,1] = gm
                coords[:, :, 1] = coords[:, :, 1] * (255 / np.max(coords[:, :, 1]))
                coords = coords.astype('uint8')
                cv2.imwrite(f'STimg/{save_name}.{i}.green.tif', coords)

        if len(args.blue) != 0:
            for i in args.blue:
                coords = np.zeros((gem_data['y'].max() + 10, gem_data['x'].max() + 10, 3), dtype=float)
                blue = gem_data[gem_data['variable']==i].copy()
                bm = blue['value'].mean() * args.contrast[2]
                coords[blue['y'], blue['x'], 0] = blue['value']
                coords[coords[:, :, 0] >= bm,0] = bm
                coords[:, :, 0] = coords[:, :, 0] * (255 / np.max(coords[:, :, 0]))
                coords = coords.astype('uint8')
                cv2.imwrite(f'STimg/{save_name}.{i}.blue.tif', coords)
        print('split plot finish')

    print('mission completed')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--function", help="set -f all/merge/split/nCount, default all", type=str, default='all')
    parser.add_argument("-R", "--red", help="-R gene1 gene2 ...", type=str, default=[], nargs='+')
    parser.add_argument("-G", "--green", help="-G gene1 gene2 ...", type=str, default=[], nargs='+')
    parser.add_argument("-B", "--blue", help="-B gene1 gene2 ...", type=str, default=[], nargs='+')
    parser.add_argument("-x", "--chop", help="x-start x-end y-start y-end", type=int, default=[], nargs='+')
    parser.add_argument("-g", "--gem", help="input the gem", type=str, default='', nargs='?')
    parser.add_argument("-b", "--bin", help="input the bin size, default 50", type=int, default=50)
    parser.add_argument("-C", "--columns", help="set the gem columns names, must include(variable,x,y,value) default: variable x y value", type=str, default=['variable', 'x', 'y', 'value'], nargs='+')
    parser.add_argument("-c", "--contrast", help="-c contrast rgb, default 3 3 3", type=int, default=[3,3,3], nargs='+')
    parser.add_argument("-n", "--normalize", help="-n 0 1000 normalize nCount, default None", type=int, default=[],nargs='+')
    parser.add_argument("-o", "--output", help="output img name, default output", type=str, default='output',nargs='?')
    args = parser.parse_args()
    main(args)



