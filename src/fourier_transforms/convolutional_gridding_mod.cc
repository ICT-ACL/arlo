#include <Python.h>
#include <ndarrayobject.h>
#include <ndarraytypes.h>
#include "convolutional_gridding.h"

static PyObject* convolutional_degrid_c(PyObject* self,PyObject* args) {

    PyArrayObject *vis_obj = NULL;
    PyArrayObject *kernels_obj = NULL;
    PyArrayObject *kernel_indices_obj = NULL;
    PyArrayObject *uvgrid_obj = NULL;
    PyArrayObject *vuvwmap_obj = NULL;
    PyArrayObject *vfrequencymap_obj = NULL;

    if (!PyArg_ParseTuple(args, "O!O!O!O!O!O!", 
                        &PyArray_Type, &vis_obj, 
                        &PyArray_Type, &kernels_obj, 
                        &PyArray_Type, &kernel_indices_obj, 
                        &PyArray_Type, &uvgrid_obj,
                        &PyArray_Type, &vuvwmap_obj,
                        &PyArray_Type, &vfrequencymap_obj)) {
        return NULL;
    }

    return Py_BuildValue("i", convolutional_degrid_c(vis_obj,
                                                    kernels_obj,
                                                    kernel_indices_obj,
                                                    uvgrid_obj,
                                                    vuvwmap_obj,
                                                    vfrequencymap_obj));
}


static PyObject* convolutional_grid_c(PyObject* self,PyObject* args) {

    PyArrayObject *uvgrid_obj = NULL;
    PyArrayObject *sumwt_obj = NULL;
    PyArrayObject *vis_obj = NULL;
    PyArrayObject *visweights_obj = NULL;
    PyArrayObject *kernels_obj = NULL;
    PyArrayObject *kernel_indices_obj = NULL;
    PyArrayObject *vuvwmap_obj = NULL;
    PyArrayObject *vfrequencymap_obj = NULL;

    if (!PyArg_ParseTuple(args, "O!O!O!O!O!O!O!O!", 
                        &PyArray_Type, &uvgrid_obj, 
                        &PyArray_Type, &sumwt_obj, 
                        &PyArray_Type, &vis_obj, 
                        &PyArray_Type, &visweights_obj, 
                        &PyArray_Type, &kernels_obj, 
                        &PyArray_Type, &kernel_indices_obj,
                        &PyArray_Type, &vuvwmap_obj,
                        &PyArray_Type, &vfrequencymap_obj)) {
        return NULL;
    }

    return Py_BuildValue("i", convolutional_grid_c(uvgrid_obj, 
                                                    sumwt_obj, 
                                                    vis_obj, 
                                                    visweights_obj, 
                                                    kernels_obj, 
                                                    kernel_indices_obj,
                                                    vuvwmap_obj,
                                                    vfrequencymap_obj));
}


static PyMethodDef module_methods[] = {
    {"convolutional_degrid_c", (PyCFunction) convolutional_degrid_c, METH_VARARGS, "calculates the fibonachi number"},
    {"convolutional_grid_c", (PyCFunction) convolutional_grid_c, METH_VARARGS, "calculates the fibonachi number"},
    {NULL, NULL, 0, NULL}
};


static struct PyModuleDef convolutional_gridding_mod = {
    PyModuleDef_HEAD_INIT,
    "convolutional_gridding_mod", /* name of module */
    "",          /* module documentation, may be NULL */
    -1,          /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    module_methods
};


PyMODINIT_FUNC PyInit_convolutional_gridding_mod(void) {
    import_array();
    return PyModule_Create(&convolutional_gridding_mod);
}