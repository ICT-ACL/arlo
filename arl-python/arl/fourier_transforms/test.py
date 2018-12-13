from fft_support_mod import fft_c, ifft_c, pad_mid_c, extract_mid_c, extract_oversampled_c
import numpy as np

if __name__ == '__main__' :

	a = np.zeros((4,4), dtype=np.complex128)
	for i in range(a.shape[0]):
		for j in range(a.shape[1]):
			a[i,j] = i + j;
	af = np.zeros(a.shape, dtype=np.complex128)
	
	fft_c(af, a)
	print(af)

	ifft_c(a, af)
	print(af)

	# pad_mid_c(af, a, 2)
	# print(af)

	# extract_mid_c(af, a, 2)
	# print(af)

	extract_oversampled_c(af, a, 8, 8, 2, 4)
	print(af)
