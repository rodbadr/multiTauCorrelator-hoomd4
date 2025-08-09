// Copyright (c) 2009-2017 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause
// License.

// Include the defined classes that are to be exported to python
// #include "Correlator.h"
#include "correlator_likh.h"
#include "correlator_pybind.h"

#include <pybind11/pybind11.h>

PYBIND11_MODULE(_multiTauCorrelator, m) {
    export_Correlator(m);
}
