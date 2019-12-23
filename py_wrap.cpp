//
// Created by Malysheva, Nadezhda on 2019-07-29.
//

#include <pybind11/pybind11.h>
#include "SSA.h"
#include "Specie.h"

using namespace std;
namespace py = pybind11;


PYBIND11_MODULE(pywrap, m)
{
    m.doc() = "Adaptive networks stochastic algorithm";
    py::class_<ContactNetwork>(m, "ContactNetwork")
            .def(py::init<size_t, size_t>())
            .def ("size", &ContactNetwork::size);
    py::class_<SSA>(m, "SSA")
            .def(py::init<>())
            .def("exe", &SSA::exe)
            .def("execute", &SSA::execute);
}
