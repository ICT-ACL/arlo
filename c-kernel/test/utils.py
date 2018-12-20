import numpy as np

def load_line(line):
    items = [int(item) for item in line.split(' ')[2:]]
    return items


def load_shape(shape_filepath):
    with open(shape_filepath) as fr:
        vshape = load_line(fr.readline())
        gshape = load_line(fr.readline())
        kshape = load_line(fr.readline())
        return vshape, gshape, kshape


def create_random_data(shape, low, high, dtype):
    if dtype == 'int':
        data = np.random.randint(low, high, shape)
    elif dtype == 'float':
        data = np.random.rand(*shape) * (high - low) + low
    elif dtype == 'complex':
        real = np.random.rand(*shape) * (high - low) + low
        imag = np.random.rand(*shape) * (high - low) + low
        data = real + imag * 1j
    return data


def store_data(filepath, arr):
    dtype = str(arr.dtype)
    if dtype.startswith('int'):
        dtype = 'int32'
    elif dtype.startswith('float'):
        dtype = 'float64'
    elif dtype.startswith('complex'):
        dtype = 'complex128'
    arr.ravel().astype(dtype).tofile(filepath)


def str2arr(string):
    return tuple(int(num) for num in string.split(','))