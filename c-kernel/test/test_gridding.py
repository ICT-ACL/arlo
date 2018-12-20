import os
import sys
sys.path.append('../../arl-python')

import numpy as np
import time
import argparse

from arl.fourier_transforms.convolutional_gridding import convolutional_grid, convolutional_degrid
from utils import *


def test_degrid(data_dir, vshape, gshape, kshape):
    kernels = []
    for i in range(kshape[0]):
        kernels.append(create_random_data(kshape[1:], 0, 1, 'complex'))
    kernel_indices = create_random_data((vshape[0]), 0, kshape[0], 'int')
    kernel_list = kernel_indices, kernels
    uvgrid = create_random_data(gshape, -1000, 1000, 'complex')
    vuvwmap = create_random_data((vshape[0], 3), -0.4, 0.4, 'float')
    vfrequencymap = create_random_data((vshape[0]), 0, 6, 'int')

    store_data(os.path.join(data_dir, 'kernel_indices.dat'), kernel_indices)
    store_data(os.path.join(data_dir, 'kernels.dat'), np.array(kernels))
    store_data(os.path.join(data_dir, 'uvgrid.dat'), uvgrid)
    store_data(os.path.join(data_dir, 'vuvwmap.dat'), vuvwmap)
    store_data(os.path.join(data_dir, 'vfrequencymap.dat'), vfrequencymap)

    start = time.time()
    vis = convolutional_degrid(kernel_list, vshape, uvgrid, vuvwmap, vfrequencymap)
    stop = time.time()
    print('Original Time:  {:.2f}s'.format(stop -start))

    store_data(os.path.join(data_dir, 'vis.dat'), vis)


def test_grid(data_dir, vshape, gshape, kshape):
    kernels = []
    for i in range(kshape[0]):
        kernels.append(create_random_data(kshape[1:], 0, 1, 'complex'))
    kernel_indices = create_random_data((vshape[0]), 0, kshape[0], 'int')
    kernel_list = kernel_indices, kernels
    uvgrid = create_random_data(gshape, -1000, 1000, 'complex')

    vis = create_random_data(vshape, -5000, 5000, 'complex')
    visweights = create_random_data(vshape, 1, 2, 'float')
    vuvwmap = create_random_data((vshape[0], 3), -0.4, 0.4, 'float')
    vfrequencymap = create_random_data((vshape[0]), 0, 6, 'int')

    store_data(os.path.join(data_dir, 'kernel_indices.dat'), kernel_indices)
    store_data(os.path.join(data_dir, 'kernels.dat'), np.array(kernels))
    store_data(os.path.join(data_dir, 'uvgrid_before.dat'), uvgrid)
    store_data(os.path.join(data_dir, 'vis.dat'), vis)
    store_data(os.path.join(data_dir, 'visweights.dat'), visweights)
    store_data(os.path.join(data_dir, 'vuvwmap.dat'), vuvwmap)
    store_data(os.path.join(data_dir, 'vfrequencymap.dat'), vfrequencymap)

    start = time.time()
    uvgrid, sumwt = convolutional_grid(kernel_list, uvgrid, vis, visweights, vuvwmap, vfrequencymap)
    stop = time.time()
    print('Original Time:  {:.2f}s'.format(stop - start))

    store_data(os.path.join(data_dir, 'uvgrid_after.dat'), uvgrid)
    store_data(os.path.join(data_dir, 'sumwt.dat'), sumwt)


if __name__ == '__main__':
    # np.random.seed(0)
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', type=str, default='grid')
    parser.add_argument('--data_dir', type=str, default='./data')
    parser.add_argument('--vshape', type=str, default='0,0')
    parser.add_argument('--gshape', type=str, default='0,0,0,0')
    parser.add_argument('--kshape', type=str, default='0,0,0,0,0')
    args = parser.parse_args()

    test = {'grid': test_grid, 'degrid': test_degrid}
    test[args.mode](args.data_dir, str2arr(args.vshape), str2arr(args.gshape), str2arr(args.kshape))