#include <Python.h>
#include <ndarrayobject.h>
#include <ndarraytypes.h>
#include "params.h"

static PyObject* get_rowmap_c(PyObject* self,PyObject* args) {

    PyArrayObject *vmap_obj = NULL;
    PyArrayObject *col_obj = NULL;
    PyArrayObject *ucol_obj = NULL;

    if (!PyArg_ParseTuple(args, "O!O!O!", 
                        &PyArray_Type, &vmap_obj, 
                        &PyArray_Type, &col_obj, 
                        &PyArray_Type, &ucol_obj)) {
        return NULL;
    }

    return Py_BuildValue("i", get_rowmap_c(vmap_obj, col_obj, ucol_obj));
}


static PyMethodDef module_methods[] = {
    {"get_rowmap_c", (PyCFunction) get_rowmap_c, METH_VARARGS, "calculates the fibonachi number"},
    {NULL, NULL, 0, NULL}
};


static struct PyModuleDef params_mod = {
    PyModuleDef_HEAD_INIT,
    "params_mod", /* name of module */
    "",          /* module documentation, may be NULL */
    -1,          /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    module_methods
};


PyMODINIT_FUNC PyInit_params_mod(void) {
    import_array();
    return PyModule_Create(&params_mod);
}