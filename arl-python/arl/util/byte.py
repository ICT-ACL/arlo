# handle byte order
import numpy as np

def native_order(x):
	if x.dtype.byteorder == '>':
		return x.byteswap()
	return x