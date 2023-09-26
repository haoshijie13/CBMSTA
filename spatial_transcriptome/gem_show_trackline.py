import os
import argparse
def main(args):
    import sys
    import time
    import numpy as np
    import cv2
    import pandas as pd
    import scipy.ndimage as nd

    def find_pattern(sums, pattern):
        pattern_len = np.sum(pattern)
        max_pattern_num = int(int((len(sums) + pattern_len - 1) // pattern_len) + 1)
        patterns = pattern * max_pattern_num
        patternarr = np.array(patterns)
        tmp_positions = np.cumsum(patternarr)
        positions = np.zeros(len(tmp_positions) + 1, dtype=int)
        positions[1:len(tmp_positions) + 1] = tmp_positions
        for shift_x in range(0, pattern_len, 1):
            tmp_indexs = positions - shift_x
            tmp_indexs = tmp_indexs[tmp_indexs >= 0]
            tmp_indexs = tmp_indexs[tmp_indexs < len(sums) - 2]
            if len(tmp_indexs) < 2:
                continue
            if np.sum(sums[tmp_indexs]) + np.sum(sums[tmp_indexs + 1]) + np.sum(sums[tmp_indexs + 2]) == 0:
                return tmp_indexs

    grid_x_715 = [112, 144, 208, 224, 224, 208, 144, 112, 160]
    grid_y_715 = [112, 144, 208, 224, 224, 208, 144, 112, 160]

    ##  chip500
    grid_x_500 = [240, 300, 330, 390, 390, 330, 300, 240, 420]
    grid_y_500 = [240, 300, 330, 390, 390, 330, 300, 240, 420]

    print(f"gem file is {args.gem}")
    print(f'chip is {args.chip}',file=sys.stderr)
    print('loading gem ...', file=sys.stderr)
    print(time.strftime("%Y-%m-%d %H:%M:%S"), file=sys.stderr, flush=True)
    gem = pd.read_csv(args.gem, sep='\s+', header=0, compression='infer', comment='#',low_memory=True,na_filter=False,engine='c')
    gem.columns = args.columns
    gem_data = gem[['variable', 'x', 'y', 'value']].copy()

    x_list = pd.unique(gem_data['x']).tolist()
    y_list = pd.unique(gem_data['y']).tolist()
    x_sum = np.zeros(max(x_list)+1)
    y_sum = np.zeros(max(y_list)+1)
    x_sum[x_list] = 1
    y_sum[y_list] = 1

    if args.chip == 'chip715':
        grid_x = grid_x_715
        grid_y = grid_y_715
    else:
        grid_x = grid_x_500
        grid_y = grid_y_500
    print("start to find trackline")
    x_indexes = find_pattern(x_sum, grid_x)
    y_indexes = find_pattern(y_sum, grid_y)
    mask = np.zeros((len(y_sum),len(x_sum)), dtype='uint8')
    mask[:, x_indexes] = 1
    mask[:, x_indexes + 1] = 1
    mask[:, x_indexes + 2] = 1
    mask[y_indexes, :] = 1
    mask[y_indexes + 1, :] = 1
    mask[y_indexes + 2, :] = 1

    print("start to combine trackline and heatmap")
    gem_data['x'] = round(gem_data['x']/args.bin, 0).astype(int)
    gem_data['y'] = round(gem_data['y']/args.bin, 0).astype(int)
    gem_data = gem_data.groupby(['variable','x', 'y']).agg(value=('value', 'sum')).reset_index()
    print("max_value: ",gem_data['value'].max()," min_value: ",gem_data['value'].min())

    expression = np.zeros((gem_data['y'].max()+1,gem_data['x'].max()+1))
    expression[gem_data['y'],gem_data['x']] = gem_data['value']

    if len(args.normalize) == 0:
        mean_contrast = expression['value'].mean() * 2
        expression[expression >= mean_contrast] = mean_contrast
        expression = expression * (255 / mean_contrast)
    else:
        expression[expression <= args.normalize[0]] = args.normalize[0]
        expression[expression >= args.normalize[1]] = args.normalize[1]
        expression = (expression-np.min(expression)) * (255 / (np.max(expression)-np.min(expression)))
    expression = expression.astype('uint8')

    scale_matrix = np.matrix([[args.bin,0,0],
                              [0,args.bin,0],
                              [0,0,1]])
    #expression = nd.affine_transform(expression.T,scale_matrix,output_shape=mask.T.shape,order=0)
    #expression= expression.T
    expression = nd.zoom(expression,args.bin,order=0)
    background = np.zeros((max(expression.shape[0],mask.shape[0]),max(expression.shape[1],mask.shape[1])))
    background[0:expression.shape[0], 0:expression.shape[1]] = expression
    background2 = np.zeros_like(background)
    background2[0:mask.shape[0], 0:mask.shape[1]] = mask
    background[background2==1]=255

    background = background.astype('uint8')
    cv2.imwrite(f'{args.output}.heatmap.marked.tif', background)
    print('mission complete')
    print(time.strftime("%Y-%m-%d %H:%M:%S"))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gem", help="gem file", type=str,default='')
    parser.add_argument("-t", "--chip", help="chip715/chip500", type=str, default='chip500')
    parser.add_argument("-b", "--bin", help="input the bin size, default 50", type=int, default=50)
    parser.add_argument("-n", "--normalize", help="-n 0 1000 normalize nCount, default None", type=int, default=[],nargs='+')
    parser.add_argument("-C", "--columns",
                        help="set the gem columns names, must include(variable,x,y,value) default: variable x y value",
                        type=str, default=['variable', 'x', 'y', 'value'], nargs='+')
    parser.add_argument("-o", "--output", help="output img name, default output", type=str, default='output')
    args = parser.parse_args()
    main(args)
