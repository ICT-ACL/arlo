#include <stdio.h>
#include "debug.h"

void print(const char *name, int *arr, int n) {
	printf("%s: [  ", name);
	if (n < 10) {
		for (int i = 0; i < n; i++) {
			printf("%d  ", arr[i]);
		}
	}
	else {
		for (int i = 0; i < 5; i++) {
			printf("%d  ", arr[i]);
		}
		printf("...  ");
		for (int i = n-5; i < n; i++) {
			printf("%d  ", arr[i]);
		}
	}
	printf("]\n");
}


void print(const char *name, double *arr, int n) {
	printf("%s: [  ", name);
	if (n < 10) {
		for (int i = 0; i < n; i++) {
			printf("%g  ", arr[i]);
		}
	}
	else {
		for (int i = 0; i < 5; i++) {
			printf("%g  ", arr[i]);
		}
		printf("...  ");
		for (int i = n-5; i < n; i++) {
			printf("%g  ", arr[i]);
		}
	}
	printf("]\n");
}


void print(const char *name, complex_t *arr, int n) {
	printf("%s: [  ", name);
	if (n < 10) {
		for (int i = 0; i < n; i++) {
			printf("%8.4lf%+8.4lfj  ", creal(arr[i]), cimag(arr[i]));
		}
	}
	else {
		for (int i = 0; i < 5; i++) {
			printf("%8.4lf%+8.4lfj  ", creal(arr[i]), cimag(arr[i]));
		}
		printf("...  ");
		for (int i = n-5; i < n; i++) {
			printf("%8.4lf%+8.4lfj  ", creal(arr[i]), cimag(arr[i]));
		}
	}
	printf("]\n");
}

void print(complex_t &x) {
	printf("%8.4lf%+8.4lfj", creal(x), cimag(x));
}

void printfull(const char *name, int *arr, int n, int split) {
	printf("%s: [  ", name);
	for (int i = 0; i < n; i++) {
		printf("%d  ", arr[i]);
		if ((i+1) % split == 0) {
			printf("\n         ");
		}
	}
	printf("]\n");
}

void printfull(const char *name, double *arr, int n, int split) {
	printf("%s: [  ", name);
	for (int i = 0; i < n; i++) {
		printf("%11.8e   ", arr[i]);
		if ((i+1) % split == 0) {
			printf("\n         ");
		}
	}
	printf("]\n");
}