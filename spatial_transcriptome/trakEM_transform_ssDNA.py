#!/usr/bin/env python3
#xxx.py file affine_matrix weight height
import numpy as np
import sys
import cv2
import scipy.ndimage as nd

def get_trackEM(param :str) -> np.matrix:
    """
        handle '-0.010963829,-0.999939895,0.999939895,-0.010963829,-129.2603788,1664.628308'

        @return reverse affine matrix.
    """
    affine = np.zeros((3,3))
    in_data = np.array(param.split(',')).astype(float).reshape((3,2))
    affine[0:2,:]=in_data.T
    affine[2] = [0,0,1]
    return np.matrix(affine).I


ssDNA = cv2.imread(sys.argv[1])
affineR = get_trackEM(sys.argv[2])
if len(ssDNA.shape) == 3:
    ssDNA = cv2.cvtColor(ssDNA, cv2.COLOR_BGR2GRAY)

affine_ssDNA = nd.affine_transform(ssDNA.T,affineR,output_shape=(int(sys.argv[3]),int(sys.argv[4])),order=0)
affine_ssDNA = affine_ssDNA.T

cv2.imwrite(f'{sys.argv[1]}.affine.tif',affine_ssDNA)




