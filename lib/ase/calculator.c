#include "calculator.h"
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>


static void create_atoms(Config *config, Input *input,
                         PyObject **pyAtomNum, PyObject **pyPos, PyObject **pyCell)
{
    int i, j, k;
    PyObject* pyTmp;
    /* atom_num */
    int count = 0;
    *pyAtomNum = PyList_New(config->tot_num);
    for (i = 0; i < config->ntype; ++i) {
        for (j = 0; j < input->nelem; ++j) {
            if (config->atom_num[i] == get_atom_num(input->atom_type[j])) {
                for (k = 0; k < config->each_num[i]; ++k) {
                    PyList_SetItem(*pyAtomNum, count, PyLong_FromLong(config->atom_num[i]));
                    count++;
                }
            }
        }
    }

    /* position */
    *pyPos = PyList_New(config->tot_num);
    for (i = 0; i < config->tot_num; ++i) {
        pyTmp = PyList_New(3);
        PyList_SetItem(pyTmp, 0, PyFloat_FromDouble(config->pos[i * 3 + 0]));
        PyList_SetItem(pyTmp, 1, PyFloat_FromDouble(config->pos[i * 3 + 1]));
        PyList_SetItem(pyTmp, 2, PyFloat_FromDouble(config->pos[i * 3 + 2]));
        PyList_SetItem(*pyPos, i, pyTmp);
    } 

    /* cell */
    *pyCell = PyList_New(3);
    for (i = 0; i < 3; ++i) {
        pyTmp = PyList_New(3);
        PyList_SetItem(pyTmp, 0, PyFloat_FromDouble(config->cell[i][0]));
        PyList_SetItem(pyTmp, 1, PyFloat_FromDouble(config->cell[i][1]));
        PyList_SetItem(pyTmp, 2, PyFloat_FromDouble(config->cell[i][2]));
        PyList_SetItem(*pyCell, i, pyTmp);
    }
}


static void ase_init(Calc *calc, Config *config, Input *input, MPI_Comm comm)
{
    if (!Py_IsInitialized()) {
        Py_Initialize();
    }

    PyObject* pySys = NULL;
    PyObject* pyPath = NULL;
    PyObject* pyName = NULL;
    PyObject* pyModule = NULL;
    PyObject* pyFunc = NULL;
    PyObject* pyArg = NULL;

    /* append path */
    pySys = PyImport_ImportModule("sys");
    pyPath = PyObject_GetAttrString(pySys, "path");
    pyName = PyUnicode_DecodeFSDefault(".");
    PyList_Append(pyPath, pyName);
    Py_XDECREF(pyName);

    char *ptr = strrchr(input->ase_calc, '/');
    char *name = (char *)malloc(sizeof(char) * 65536);
    if (ptr != NULL) {
        char path[65536];
        strcpy(path, input->ase_calc);
        path[strlen(input->ase_calc) - strlen(ptr)] = '\0';
        pyName = PyUnicode_DecodeFSDefault(path);
        PyList_Append(pyPath, pyName);
        Py_XDECREF(pyName);
        strcpy(name, ptr + 1);
    } else {
        strcpy(name, input->ase_calc); 
    }

    /* import */
    name[strlen(name) - 3] = '\0';
    pyName = PyUnicode_DecodeFSDefault(name);
    free(name);
    if (pyName != NULL) {
        pyModule = PyImport_Import(pyName);
        Py_XDECREF(pyName);
    } else {
        if (PyErr_Occurred()) {
            PyErr_Print();
        }
    }

    if (pyModule != NULL) {
        pyFunc = PyObject_GetAttrString(pyModule, "ase_initialize");
        if ((pyFunc != NULL) && (PyCallable_Check(pyFunc))) {
            pyArg = PyTuple_New(2);
            PyTuple_SetItem(pyArg, 0, PyUnicode_FromString(input->model_path));
            if (input->device != NULL) {
                int rank;
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                char *device = (char *)malloc(sizeof(char) * 16);
                strcpy(device, input->device);
                char *index = (char *)malloc(sizeof(char) * 8);
                sprintf(index, ":%d", rank);
                strcat(device, index);
                PyTuple_SetItem(pyArg, 1, PyUnicode_FromString(device));
                free(device);
                free(index);
            } else {
                PyTuple_SetItem(pyArg, 1, NULL);
            }
            PyObject_CallObject(pyFunc, pyArg);
            if (PyErr_Occurred()) {
                PyErr_Print();
            }
            Py_XDECREF(pyArg);
            Py_XDECREF(pyFunc);
        } else {
            if (PyErr_Occurred()) {
                PyErr_Print();
            }
        }
    } else {
        if (PyErr_Occurred()) {
            PyErr_Print();
        }
    }

    Py_XDECREF(pySys);
    Py_XDECREF(pyPath);

    /* turn on initialized tag */
    calc->initialized = 1;
    calc->ase = (void *)pyModule;
}


static void run_ase(Calc *calc, char *name, PyObject *pyArg,
                    double *output_scalar, double *output_array)
{
    int i, j;
    PyObject* pyFunc = NULL;
    PyObject* pyValue = NULL;
    PyObject* pyScalar = NULL;
    PyObject* pyVector = NULL;
    PyObject* pyMatrix = NULL;

    /* run */
    pyFunc = PyObject_GetAttrString((PyObject *)calc->ase, name);
    if (pyFunc != NULL && PyCallable_Check(pyFunc)) {
        pyValue = PyObject_CallObject(pyFunc, pyArg);
        Py_XDECREF(pyArg);
        if (pyValue != NULL && PyTuple_Check(pyValue)) {
            pyScalar = PyTuple_GetItem(pyValue, 0);
            if (pyScalar != NULL && PyFloat_Check(pyScalar)) {
                *output_scalar = PyFloat_AsDouble(pyScalar);
            } else {
                if (PyErr_Occurred()) {
                    PyErr_Print();
                }
            }
            pyMatrix = PyTuple_GetItem(pyValue, 1);
            if (pyMatrix != NULL && PyList_Check(pyMatrix)) {
                for (i = 0; i < PyList_Size(pyMatrix); ++i) {
                    pyVector = PyList_GetItem(pyMatrix, i);
                    if (pyVector != NULL && PyList_Check(pyVector)) {
                        for (j = 0; j < 3; ++j) {
                            pyScalar = PyList_GetItem(pyVector, j);
                            if (pyScalar != NULL && PyFloat_Check(pyScalar)) {
                                output_array[i * 3 + j] = PyFloat_AsDouble(pyScalar);
                            } else {
                                if (PyErr_Occurred()) {
                                    PyErr_Print();
                                }
                            }
                        }
                    } else {
                        if (PyErr_Occurred()) {
                            PyErr_Print();
                        }
                    }
                } 
            } else {
                if (PyErr_Occurred()) {
                    PyErr_Print();
                }
            }
            Py_XDECREF(pyValue);
        } else {
            if (PyErr_Occurred()) {
                PyErr_Print();
            }
        }
        Py_XDECREF(pyFunc);
    } else {
        if (PyErr_Occurred()) {
            PyErr_Print();
        }
    }
}


void oneshot(Calc *calc, Config *config, Input *input,
             double *energy, double *force, MPI_Comm comm)
{
    PyObject* pyArg = NULL;
    PyObject* pyAtomNum = NULL;
    PyObject* pyPos = NULL;
    PyObject* pyCell = NULL;

    /* initialize */
    if (calc->initialized == 0) {
        ase_init(calc, config, input, comm);
    }
    create_atoms(config, input, &pyAtomNum, &pyPos, &pyCell);

    /* argument */
    pyArg = PyTuple_New(3);
    PyTuple_SetItem(pyArg, 0, pyAtomNum);
    PyTuple_SetItem(pyArg, 1, pyPos);
    PyTuple_SetItem(pyArg, 2, pyCell);

    run_ase(calc, "oneshot", pyArg, energy, force);
}


void atom_relax(Calc *calc, Config *config, Input *input,
                double *energy, MPI_Comm comm)
{
    int i;
    PyObject* pyArg;
    PyObject* pyAtomNum;
    PyObject* pyPos;
    PyObject* pyCell;
    PyObject* pyFix;

    /* initialize */
    if (calc->initialized == 0) {
        ase_init(calc, config, input, comm);
    }
    create_atoms(config, input, &pyAtomNum, &pyPos, &pyCell);

    /* fix */
    pyFix = PyList_New(config->tot_num);
    for (i = 0; i < config->tot_num; ++i) {
        if (config->fix[i] == 1) {
            PyList_SetItem(pyFix, i, Py_True);
        } else {
            PyList_SetItem(pyFix, i, Py_False);
        }
    }

    /* argument */
    pyArg = PyTuple_New(5);
    PyTuple_SetItem(pyArg, 0, pyAtomNum);
    PyTuple_SetItem(pyArg, 1, pyPos);
    PyTuple_SetItem(pyArg, 2, pyCell);
    PyTuple_SetItem(pyArg, 3, pyFix);
    PyTuple_SetItem(pyArg, 4, PyFloat_FromDouble(input->f_tol));

    run_ase(calc, "atom_relax", pyArg, energy, config->pos);
}


void free_calc(Calc *calc)
{
    Py_XDECREF((PyObject *)calc->ase);
    /* delete python instance */
    if (Py_IsInitialized()) {
        Py_Finalize();
    }
}
