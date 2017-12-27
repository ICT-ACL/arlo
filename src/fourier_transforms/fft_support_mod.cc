#include <Python.h>
#include <ndarrayobject.h>
#include <ndarraytypes.h>
#include "fft_support.h"

static PyObject* fft_c(PyObject* self,PyObject* args) {

    PyArrayObject *af_obj = NULL;
    PyArrayObject *a_obj = NULL;

    if (!PyArg_ParseTuple(args, "O!O!", 
                        &PyArray_Type, &af_obj, 
                        &PyArray_Type, &a_obj)) {
        return NULL;
    }

    return Py_BuildValue("i", fft_c(af_obj, a_obj));
}


static PyObject* ifft_c(PyObject* self,PyObject* args) {

    PyArrayObject *a_obj = NULL;
    PyArrayObject *af_obj = NULL;

    if (!PyArg_ParseTuple(args, "O!O!", 
                        &PyArray_Type, &a_obj, 
                        &PyArray_Type, &af_obj)) {
        return NULL;
    }

    return Py_BuildValue("i", ifft_c(a_obj, af_obj));
}


static PyObject* pad_mid_c(PyObject* self,PyObject* args) {

    PyArrayObject *out_obj = NULL;
    PyArrayObject *in_obj = NULL;
    int npixel;

    if (!PyArg_ParseTuple(args, "O!O!i", 
                        &PyArray_Type, &out_obj, 
                        &PyArray_Type, &in_obj,
                        &npixel)) {
        return NULL;
    }

    return Py_BuildValue("i", pad_mid_c(out_obj, in_obj, npixel));
}


static PyObject* extract_mid_c(PyObject* self,PyObject* args) {

    PyArrayObject *out_obj = NULL;
    PyArrayObject *in_obj = NULL;
    int npixel;

    if (!PyArg_ParseTuple(args, "O!O!i", 
                        &PyArray_Type, &out_obj, 
                        &PyArray_Type, &in_obj,
                        &npixel)) {
        return NULL;
    }

    return Py_BuildValue("i", extract_mid_c(out_obj, in_obj, npixel));
}


static PyObject* extract_oversampled_c(PyObject* self,PyObject* args) {

    PyArrayObject *mid_obj = NULL;
    PyArrayObject *a_obj = NULL;
    int xf, yf;
    int kernel_oversampling;
    int kernelwidth;

    if (!PyArg_ParseTuple(args, "O!O!iiii", 
                        &PyArray_Type, &mid_obj, 
                        &PyArray_Type, &a_obj,
                        &xf,
                        &yf,
                        &kernel_oversampling,
                        &kernelwidth)) {
        return NULL;
    }

    return Py_BuildValue("i", extract_oversampled_c(mid_obj, a_obj, xf, yf, 
                                                    kernel_oversampling, kernelwidth));
}


static PyMethodDef module_methods[] = {
    {"fft_c", (PyCFunction) fft_c, METH_VARARGS, "calculates the fibonachi number"},
    {"ifft_c", (PyCFunction) ifft_c, METH_VARARGS, "calculates the fibonachi number"},
    {"pad_mid_c", (PyCFunction) pad_mid_c, METH_VARARGS, "calculates the fibonachi number"},
    {"extract_mid_c", (PyCFunction) extract_mid_c, METH_VARARGS, "calculates the fibonachi number"},
    {"extract_oversampled_c", (PyCFunction) extract_oversampled_c, METH_VARARGS, "calculates the fibonachi number"},
    {NULL, NULL, 0, NULL}
};


static struct PyModuleDef fft_support_mod = {
    PyModuleDef_HEAD_INIT,
    "fft_support_mod", /* name of module */
    "",          /* module documentation, may be NULL */
    -1,          /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    module_methods
};


PyMODINIT_FUNC PyInit_fft_support_mod(void) {
    import_array();
    return PyModule_Create(&fft_support_mod);
}