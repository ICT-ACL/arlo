#include <map>
#include <set>
#include <stdlib.h>
#include <math.h>

#include "params.h"

using namespace std;

void get_rowmap(int *vmap, double *col, double *ucol, int len, int ulen) {

    if (ucol == NULL) {
        set<double> colset;
        for (int i = 0; i < len; i++) {
            colset.insert(col[i]);
        }

        ulen = int(colset.size());
        ucol = (double*)malloc(sizeof(double) * ulen);
        int i = 0;
        for (set<double>::iterator it = colset.begin(); it != colset.end(); ++it) {
            ucol[i++] = *it;
        }
    }

    map<int, int> pdict;
    for (int i = 0; i < ulen; i++) {
        int val = int(round(ucol[i]));
        pdict[val] = i;
    }

    for (int i = 0; i < len; i++) {
        int val = int(round(col[i]));
        vmap[i] = pdict[val];
    }
}


int get_rowmap_c(PyArrayObject *&vmap_obj,
                 PyArrayObject *&col_obj,
                 PyArrayObject *&ucol_obj) {

    int len = col_obj->dimensions[0];
    int ulen = ucol_obj->dimensions[0];

    int *vmap = (int*) vmap_obj->data;
    double *col = (double*)col_obj->data;
    double *ucol = (double*)ucol_obj->data;

    get_rowmap(vmap, col, ucol, len, ulen);

    return 0;
}