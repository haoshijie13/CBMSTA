#!/usr/bin/env python3
import argparse
import os

def main(args):
    import pandas as pd
    from skimage import io
    print('script start')
    gem_data = pd.read_csv(args.gem, sep='\s+', header=0, compression='infer', comment='#',low_memory=True,na_filter=False,engine='c')
    gem_data['bx'] = round(gem_data['x'] / args.bin, 0).astype(int)
    gem_data['by'] = round(gem_data['y'] / args.bin, 0).astype(int)
    if args.mask != '':
        mask = io.imread(args.mask)
        gem_data["mask"] = mask[gem_data['by'],gem_data['bx']]
        gem_data = gem_data[gem_data['mask']!=0]
    if len(args.label[0]) != 0:
        gem_data['bcoor'] = gem_data['bx'].map(str) + '_' + gem_data['by'].map(str)
        label = pd.read_csv(args.label[0],header=0,sep='\t',low_memory=False)
        label['coor'] = label[args.label[1]].map(str) + '_' + label[args.label[2]].map(str)
        gem_data = gem_data[gem_data['bcoor'].isin(label['coor'])]
        for i in args.label[3:]:
            gem_data[i] = gem_data['bcoor'].map(dict(zip(label['coor'],label[i])))
        gem_data.drop(columns='bcoor', inplace=True)

    gem_data.drop(columns='bx',inplace = True)
    gem_data.drop(columns='by',inplace = True)
    if not os.path.exists("gem_cut"):
        print("create gem_cut")
        os.mkdir("gem_cut")
    gem_data.to_csv(f'gem_cut/{args.output}.cut.gem.gz', sep='\t', header=True, index=False)
    print(f'file save in gem_cut/{args.output}.cut.gem.gz')
    print('mission completed')




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gem", help="input the gem file path", type=str, default='')
    parser.add_argument("-m", "--mask", help="input the mask file path", type=str, default='')
    parser.add_argument("-l", "--label", help="-l filepath x(column) y(column) column1 column2;input the label file path and label columns name", type=str, nargs='+', default=[])
    parser.add_argument("-b", "--bin", help="input the bin size, default 50", type=int, default=50)
    parser.add_argument("-o", "--output", help="output img name, default output", type=str, default='output',nargs='?')
    args = parser.parse_args()
    main(args)
