// image/cleaners_mod.cc
// cleaners python module
// Author: You Haihang, Yang Runkai, Liu Tao,
// Applied Computing Lab, ICT, CAS 

#include <Python.h>
#include <ndarrayobject.h>
#include <ndarraytypes.h>
#include "cleaners.h"

static PyObject* msmfsclean_kernel_c(PyObject* self,PyObject* args) {

    PyArrayObject *m_model_obj = NULL;
    PyArrayObject *residual_obj = NULL;
    PyArrayObject *scalestack_obj = NULL;
    PyArrayObject *smresidual_obj = NULL;
    PyArrayObject *ssmmpsf_obj = NULL;
    PyArrayObject *hsmmpsf_obj = NULL;
    PyArrayObject *ihsmmpsf_obj = NULL;
    PyArrayObject *ldirty_obj = NULL;
    PyArrayObject *psf_obj = NULL;
    PyArrayObject *windowstack_obj = NULL;

    double gain;
    double absolutethresh;
    int niter;
    char* findpeak;

    if (!PyArg_ParseTuple(args, "O!O!O!O!O!O!O!O!O!O!ddis", 
                        &PyArray_Type, &m_model_obj,
                        &PyArray_Type, &residual_obj,
                        &PyArray_Type, &scalestack_obj,
                        &PyArray_Type, &smresidual_obj,
                        &PyArray_Type, &ssmmpsf_obj,
                        &PyArray_Type, &hsmmpsf_obj,
                        &PyArray_Type, &ihsmmpsf_obj,
                        &PyArray_Type, &ldirty_obj,
                        &PyArray_Type, &psf_obj,
                        &PyArray_Type, &windowstack_obj,
                        &gain, &absolutethresh, &niter, &findpeak)) {
        return NULL;
    }

    return Py_BuildValue("i", msmfsclean_kernel_c(m_model_obj,
                                                residual_obj,
                                                scalestack_obj,
                                                smresidual_obj,
                                                ssmmpsf_obj,
                                                hsmmpsf_obj,
                                                ihsmmpsf_obj,
                                                ldirty_obj,
                                                psf_obj,
                                                windowstack_obj,
                                                gain,
                                                absolutethresh,
                                                niter,
                                                findpeak));
}


static PyMethodDef module_methods[] = {
    {"msmfsclean_kernel_c", (PyCFunction) msmfsclean_kernel_c, METH_VARARGS, "calculates the fibonachi number"},
    {NULL, NULL, 0, NULL}
};


static struct PyModuleDef cleaners_mod = {
    PyModuleDef_HEAD_INIT,
    "cleaners_mod", /* name of module */
    "",          /* module documentation, may be NULL */
    -1,          /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    module_methods
};


PyMODINIT_FUNC PyInit_cleaners_mod(void) {
    import_array();
    return PyModule_Create(&cleaners_mod);
}
