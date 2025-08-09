#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "correlator_likh.h"
#include "correlator_pybind.h"

void export_Correlator(pybind11::module& m) {
  pybind11::class_<Correlator_Likh>(m, "Correlator_Likh")
        .def(pybind11::init<int, int, int>())
        .def("initialize", &Correlator_Likh::initialize)
        .def("add", &Correlator_Likh::add, pybind11::arg("w"), pybind11::arg("k") = 0)
        .def("evaluate", &Correlator_Likh::evaluate, pybind11::arg("norm") = false)
        .def_property_readonly("npcorr", [](const Correlator_Likh& self) { return self.npcorr; })
        .def_property_readonly("t", [](const Correlator_Likh& self) {
            return pybind11::array_t<double>(self.npcorr, self.t);
        })
        .def_property_readonly("f", [](const Correlator_Likh& self) {
            return pybind11::array_t<double>(self.npcorr, self.f);
        })
        ;
}
