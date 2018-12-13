import numpy as np

def add_arr_sum(*arrs):
    add_sum_arr = np.array(arrs)
    add_sum_arr = np.around(np.append(add_sum_arr, add_sum_arr.sum()), decimals=2)
    return add_sum_arr