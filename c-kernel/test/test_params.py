import os
import sys

sys.path.append('../../arl-python')

import time
import numpy as np
import argparse

from arl.imaging.params import get_rowmap
from utils import store_data

def test_params(data_dir, nrows, low, high, len):
    frequency = np.random.choice(np.linspace(low, high, len), nrows)
    ufrequency = np.unique(frequency)

    store_data(os.path.join(data_dir, 'frequency.dat'), frequency)
    store_data(os.path.join(data_dir, 'ufrequency.dat'), ufrequency)

    start = time.time()
    vmap = np.array(get_rowmap(frequency, ufrequency))
    stop = time.time()
    print('Original Time:  {:.2f}s'.format(stop - start))

    store_data(os.path.join(data_dir, 'vmap.dat'), vmap)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_dir', type=str, default='./data')
    parser.add_argument('--nrows', type=int, default=0)
    parser.add_argument('--low', type=float, default=0)
    parser.add_argument('--high', type=float, default=0)
    parser.add_argument('--len', type=int, default=0)
    args = parser.parse_args()

    test_params(args.data_dir, args.nrows, args.low, args.high, args.len)
