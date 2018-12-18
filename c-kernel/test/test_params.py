import os
import sys

sys.path.append('../../arl-python')

import time
import numpy as np

from arl.imaging.params import get_rowmap
from utils import store_data

def test_params(data_dir, nrows):
    frequency = np.random.choice(np.linspace(0.8e8, 1.2e8, 7), nrows)
    ufrequency = np.unique(frequency)

    store_data(os.path.join(data_dir, 'frequency.dat'), frequency)
    store_data(os.path.join(data_dir, 'ufrequency.dat'), ufrequency)

    start = time.time()
    vmap = np.array(get_rowmap(frequency, ufrequency))
    stop = time.time()

    store_data(os.path.join(data_dir, 'vmap.dat'), vmap)

    print('Original Time:  {:.2f}s'.format(stop - start))


if __name__ == '__main__':
    data_dir = './data'
    test_params(data_dir, nrows=479325)