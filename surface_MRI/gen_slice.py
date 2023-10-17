import os
import time
import numpy as np
import pandas as pd 
import nibabel as nib
from numba import jit
from tqdm import tqdm


def scale_xy(xy_raw, scale):
    return np.round(xy_raw / scale).astype(int)

@jit(nopython = True)
def totalRNA(xy, count, mat):
    for coor, c in zip(xy, count):
        mat[coor[0], coor[1]] += c
    return mat



def time_elapsed(start, str):
    print('{}: {:.2f}s'.format(str, time.time() - start))
    return time.time()

def gen_slice(path_input, path_output, scale_factor=50):

    start = time.time()
    print('{}:'.format(path_input))
    data_raw_txt = pd.read_csv(path_input, delimiter='\t', comment='#', header=0)
    last_end = time_elapsed(start, 'load_txt')

    MIDCount = np.array(data_raw_txt.iloc[:, 3])
    xy_raw = np.array(data_raw_txt.iloc[:, 1:3])
    xy_scaled = scale_xy(xy_raw, scale_factor)
    x_max, y_max = xy_scaled[:, 0].max() + 1, xy_scaled[:, 1].max() + 1
    MIDCount_Matrix = np.zeros((x_max, y_max))
    Mat = totalRNA(xy_scaled, MIDCount, MIDCount_Matrix)
    last_end = time_elapsed(last_end, 'gen_mat')

    nifti_img = nib.Nifti1Image(Mat, np.eye(4))
    nib.save(nifti_img, path_output)
    time_elapsed(last_end, 'save_nifti')
    print('\n')
    return True


def rescale(path_nii_slice, path_output, edge=1000):
    raw = nib.load(path_nii_slice)
    slice = raw.get_fdata()
    slice_shape = slice.shape
    new_slice = np.zeros((edge, edge))
    new_slice[:slice_shape[0], :slice_shape[1]] = slice[:1000, :1000]
    img = nib.Nifti1Image(new_slice, np.eye(4))
    nib.save(img, path_output)


def nii_3Dto2D(path_nii_volume, path_output, slice_dim=2, prefix=''):
    raw = nib.load(path_nii_volume)
    volume = raw.get_fdata()
    shape = list(range(volume.ndim))
    slice_shape = list(volume.shape)
    slice_shape[slice_dim] = 1
    new_shape = [shape.pop(slice_dim)] + shape
    volume = volume.transpose(new_shape)
    for i, slice in tqdm(enumerate(volume)):
        slice.reshape(slice_shape)
        img = nib.Nifti1Image(slice, np.eye(4))
        path_slice = os.path.join(path_output, '{}_{}.nii.gz'.format(prefix, str(i).zfill(2)))
        nib.save(img, path_slice)


def flip_rotate(path_nii_slice, path_output, manipulation=0):
    raw = nib.load(path_nii_slice)
    slice = raw.get_fdata()
    if manipulation != 0:
        if manipulation < 0:
            slice = np.filp(slice, axis=1)
        slice = np.rot90(slice, abs(manipulation))
    img = nib.Nifti1Image(slice, np.eye(4))
    nib.save(img, path_output)


@jit(nopython = True)
def xy_get_mask(xy, mask, xy_new_property):
    for i, coor in enumerate(xy):
        if coor[0] < 1600 and coor[1] < 1600:
            xy_new_property[i] = mask[coor[0], coor[1]]
    return xy_new_property

def gen_txt(path_raw, path_mask, path_output, new_property='mask', scale_factor=50):
    
    start = time.time()
    print('{}:'.format(path_raw))
    data = pd.read_csv(path_raw, delimiter='\t', comment='#', header=0)
    last_end = time_elapsed(start, 'load_txt')

    mask_full = nib.load(path_mask)
    mask = mask_full.get_fdata()
    last_end = time_elapsed(last_end, 'load_nii')

    xy_raw = np.array(data.iloc[:, 1:3])
    xy_scaled = scale_xy(xy_raw, scale_factor)
    xy_new_property = np.zeros(xy_scaled.shape[0], dtype=int)
    xy_new_property = xy_get_mask(xy_scaled, np.squeeze(mask.astype(int)), xy_new_property)
    time_elapsed(last_end, 'map_new_property')

    data[new_property] = xy_new_property
    data.to_csv(path_output, sep='\t', header=True, index=False)
    time_elapsed(last_end, 'save_txt')
    print('\n')
    return True


def gen_slice_masked(path_input, path_output, scale_factor=50):

    start = time.time()
    print('{}:'.format(path_input))
    data_raw_txt = pd.read_csv(path_input, delimiter='\t', comment='#', header=0)
    last_end = time_elapsed(start, 'load_txt')

    MIDCount = np.array(data_raw_txt.iloc[:, 3]) * np.array(data_raw_txt.iloc[:, 4])
    xy_raw = np.array(data_raw_txt.iloc[:, 1:3])
    xy_scaled = scale_xy(xy_raw, scale_factor)
    x_max, y_max = xy_scaled[:, 0].max() + 1, xy_scaled[:, 1].max() + 1
    MIDCount_Matrix = np.zeros((x_max, y_max))
    Mat = totalRNA(xy_scaled, MIDCount, MIDCount_Matrix)
    last_end = time_elapsed(last_end, 'gen_mat')

    nifti_img = nib.Nifti1Image(Mat, np.eye(4))
    nib.save(nifti_img, path_output)
    time_elapsed(last_end, 'save_nifti')
    print('\n')
    return True

def make_ROI(xy, count, mat):
    for coor, c in zip(xy, count):
        #if coor[1] > -5000 :
        mat[coor[0], coor[1]] = c
    return mat


def make_ROI_repeat(xy, count, mat):
    for coor, c in zip(xy, count):
        #if coor[1] > -5000 :
        mat[coor[0], coor[1]] = mat[coor[0], coor[1]] + c
    return mat


def gen_slice_ROI(path_input, path_output, scale_factor=50):

    start = time.time()
    print('{}:'.format(path_input))
    data_raw_txt = pd.read_csv(path_input, delimiter='\t', comment='#', header=0)
    last_end = time_elapsed(start, 'load_txt')

    MIDCount = np.array(data_raw_txt.iloc[:, 4]) #the col of the ROI in txt
    xy_raw = np.array(data_raw_txt.iloc[:, 1:3]) # the col of the xy in txt
    xy_scaled = scale_xy(xy_raw, scale_factor)
    x_max, y_max = xy_scaled[:, 0].max() + 1, xy_scaled[:, 1].max() + 1
    MIDCount_Matrix = np.zeros((x_max, y_max))

    Mat = make_ROI(xy_scaled, MIDCount, MIDCount_Matrix)
    last_end = time_elapsed(last_end, 'gen_mat')

    nifti_img = nib.Nifti1Image(Mat, np.eye(4))
    nib.save(nifti_img, path_output)
    time_elapsed(last_end, 'save_nifti')
    print('\n')
    return True


def gen_slice_ROI_from_parquet(path_input, path_output, scale_factor=50,cell_type='Purkinje'):

    start = time.time()
    print('{}:'.format(path_input))
    data_raw_txt = pd.read_parquet(path_input)
    last_end = time_elapsed(start, 'load_txt')

    data_raw_txt['cell_type'] =  data_raw_txt['cell_type'].fillna(0)
    data_raw_txt.loc[ data_raw_txt['cell_type'] == cell_type, 'cell_type'] = 1

    MIDCount = np.array(data_raw_txt.iloc[:, 5]) #the col of the ROI in txt
    xy_raw = np.array(data_raw_txt.iloc[:, 1:3]) # the col of the xy in txt
    xy_scaled = scale_xy(xy_raw, scale_factor)
    x_max, y_max = xy_scaled[:, 0].max() + 1, xy_scaled[:, 1].max() + 1
    MIDCount_Matrix = np.zeros((x_max, y_max))

    Mat = make_ROI_repeat(xy_scaled, MIDCount, MIDCount_Matrix)
    last_end = time_elapsed(last_end, 'gen_mat')

    nifti_img = nib.Nifti1Image(Mat, np.eye(4))
    nib.save(nifti_img, path_output)
    time_elapsed(last_end, 'save_nifti')
    print('\n')
    return True





def gen_rxy_slice_ROI(path_input, path_output, scale_factor=50):

    start = time.time()
    print('{}:'.format(path_input))
    data_raw_txt = pd.read_csv(path_input, delimiter='\t', comment='#', header=0)
    last_end = time_elapsed(start, 'load_txt')

    MIDCount = np.array(data_raw_txt.iloc[:, 5]) #the col of the ROI in txt
    xy_raw = np.array(data_raw_txt.iloc[:, 6:8]) # the col of the rxy in txt
    xy_scaled = scale_xy(xy_raw, scale_factor)
    x_max, y_max = xy_scaled[:, 0].max() + 1, xy_scaled[:, 1].max() + 1
    MIDCount_Matrix = np.zeros((x_max, y_max))

    Mat = make_ROI(xy_scaled, MIDCount, MIDCount_Matrix)
    last_end = time_elapsed(last_end, 'gen_mat')

    nifti_img = nib.Nifti1Image(Mat, np.eye(4))
    nib.save(nifti_img, path_output)
    time_elapsed(last_end, 'save_nifti')
    print('\n')
    return True


def gen_rxy_slice(path_input, path_output, scale_factor=50):

    start = time.time()
    print('{}:'.format(path_input))
    data_raw_txt = pd.read_csv(path_input, delimiter='\t', comment='#', header=0)
    last_end = time_elapsed(start, 'load_txt')

    MIDCount = np.array(data_raw_txt.iloc[:, 3])#the col of the midcount in txt
    xy_raw = np.array(data_raw_txt.iloc[:, 6:8])# the col of the rxy in txt
    xy_scaled = scale_xy(xy_raw, scale_factor)
    x_max, y_max = xy_scaled[:, 0].max() + 1, xy_scaled[:, 1].max() + 1
    MIDCount_Matrix = np.zeros((x_max, y_max))
    Mat = totalRNA(xy_scaled, MIDCount, MIDCount_Matrix)
    last_end = time_elapsed(last_end, 'gen_mat')

    nifti_img = nib.Nifti1Image(Mat, np.eye(4))
    nib.save(nifti_img, path_output)
    time_elapsed(last_end, 'save_nifti')
    print('\n')
    return True



