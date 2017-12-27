#ifndef PARAMS_H
#define PARAMS_H

#include <ndarrayobject.h>

void get_rowmap(int *vmap, double *col, double *ucol, int len, int ulen);

int get_rowmap_c(PyArrayObject *&vmap_obj, 
				 PyArrayObject *&col_obj,
				 PyArrayObject *&ucol_obj);

#endif