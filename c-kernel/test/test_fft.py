import os
import sys
sys.path.append('../../arl-python')

import numpy as np
import time

from arl.fourier_transforms.fft_support import fft, ifft
from utils import *


def test_fft(data_dir):
    a = create_random_data((7, 2, 1024, 1024), -100, 100, 'complex')

    start = time.time()
    ia = fft(a)
    stop = time.time()
    print('Original Time:  {:.2f}s'.format(stop - start))

    store_data(os.path.join(data_dir, 'a.dat'), a)
    store_data(os.path.join(data_dir, 'ia.dat'), ia)


def test_ifft(data_dir):
    ia = create_random_data((7, 2, 1024, 1024), -100, 100, 'complex')

    start = time.time()
    a = ifft(ia)
    stop = time.time()
    print('Original Time:  {:.2f}s'.format(stop - start))

    store_data(os.path.join(data_dir, 'ia.dat'), ia)
    store_data(os.path.join(data_dir, 'a.dat'), a)


if __name__ == '__main__':
    # np.random.seed(0)
    data_dir = './data/'

    if sys.argv[1] == 'fft':
        test_fft(data_dir)
    elif sys.argv[1] == 'ifft':
        test_ifft(data_dir)