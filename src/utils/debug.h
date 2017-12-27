# ifndef DEBUG_H
# define DEBUG_H

#include "fftw.h"

void print(const char *name, int *arr, int n);
void print(const char *name, double *arr, int n);
void print(const char *name, complex_t *arr, int n);

void print(complex_t &x);

void printfull(const char *name, int *arr, int n, int split);
void printfull(const char *name, double *arr, int n, int split);

# endif