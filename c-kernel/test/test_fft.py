import os
import sys
sys.path.append('../../arl-python')

import numpy as np
import time
import argparse

from arl.fourier_transforms.fft_support import fft, ifft
from utils import *


def test_fft(data_dir, shape):
    a = create_random_data(shape, -100, 100, 'complex')

    start = time.time()
    ia = fft(a)
    stop = time.time()
    print('Original Time:  {:.2f}s'.format(stop - start))

    store_data(os.path.join(data_dir, 'a.dat'), a)
    store_data(os.path.join(data_dir, 'ia.dat'), ia)


def test_ifft(data_dir, shape):
    ia = create_random_data(shape, -100, 100, 'complex')

    start = time.time()
    a = ifft(ia)
    stop = time.time()
    print('Original Time:  {:.2f}s'.format(stop - start))

    store_data(os.path.join(data_dir, 'ia.dat'), ia)
    store_data(os.path.join(data_dir, 'a.dat'), a)


if __name__ == '__main__':
    # np.random.seed(0)

    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', type=str, default='fft')
    parser.add_argument('--data_dir', type=str, default='./data')
    parser.add_argument('--shape', type=str, default='0,0,0,0')
    args = parser.parse_args()

    test = {'fft': test_fft, 'ifft': test_ifft}
    test[args.mode](args.data_dir, str2arr(args.shape))