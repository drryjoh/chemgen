#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "chem.hpp"

namespace py = pybind11;

PYBIND11_MODULE(chemwrapper, m) {
    m.def("source_jacobian", &source_jacobain, "Scale species vector by a factor");
}

