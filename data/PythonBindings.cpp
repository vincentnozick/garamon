// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// Copyright (c) 2018 by Norwegian University of Science and Technology
// PythonBindings.cpp
// This file is part of the Garamon for project_namespace.
// Authors: Stephane Breuils, Vincent Nozick, and Lars Tingelstad
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

/// \file PythonBindings.cpp
/// \author Lars Tingelstad
/// \brief Python bindings using pybind11.

#include "project_namespace/Mvec.hpp"

#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/iostream.h>

namespace py = pybind11;


/*!
 * @namespace project_namespace
 */
namespace project_namespace {

PYBIND11_MODULE(project_namespace_py, m) {

  project_basis_vector_index
  
  m.def("metric", [](){return metric;});

  m.attr("scalar") = 0;
  project_static_multivector_one_component_python

  // Class definition
  auto mvec = py::class_<Mvec<double>>(m, "Mvec");
  // Constructors
  mvec.def(py::init<>());
  // Get/Set
  mvec.def("__setitem__",
           [](Mvec<double>& mv, int idx, double value) { mv[idx] = value; });
  mvec.def("__getitem__",
           [](Mvec<double>& mv, int idx) { return mv[idx]; });
  // Operators
  mvec.def(py::self + py::self)
      .def(py::self + float())
      .def(float() + py::self)
      .def(py::self += py::self)
      .def(py::self - py::self)
      .def(py::self - float())
      .def(float() - py::self)
      .def(py::self -= py::self)
      .def(py::self * py::self)
      .def(py::self * float())
      .def(float() * py::self)
      .def(py::self *= py::self)
      .def(py::self / py::self)
      .def(py::self / float())
      .def(float() / py::self)
      .def(py::self /= py::self)
      .def(py::self < py::self)
      .def(py::self < float())
      .def(float() < py::self)
      .def(py::self > py::self)
      .def(py::self > float())
      .def(float() > py::self)
      .def("__invert__",
           [](const Mvec<double>& a) { return ~a; })
      .def("__eq__",
           [](Mvec<double>& a, const Mvec<double>& b) { return a == b; })
      .def("__neq__",
           [](Mvec<double>& a, const Mvec<double>& b) { return a != b; })
      .def("__or__",
           [](const Mvec<double>& a, const Mvec<double>& b) { return a | b; })
      .def("__or__",
           [](const Mvec<double>& a, double b) { return a ^ b; })
      .def("__ror__",
           [](const Mvec<double>& a, double b) { return a ^ b; })
      .def("__ior__",
           [](Mvec<double>& a, const Mvec<double>& b) { a ^= b; return a; })
      .def("__xor__",
           [](const Mvec<double>& a, const Mvec<double>& b) { return a ^ b; })
      .def("__xor__",
           [](const Mvec<double>& a, double b) { return a ^ b; })
      .def("__rxor__",
           [](const Mvec<double>& a, double b) { return a ^ b; })
      .def("__ixor__",
           [](Mvec<double>& a, const Mvec<double>& b) { a ^= b; return a; });

  // Print
  mvec.def("__repr__",
           [](const Mvec<double>& mv) {
             std::stringstream ss;
             ss << mv;
             return ss.str();
           })
      .def("norm", &Mvec<double>::norm)
      .def("quadratic_norm", &Mvec<double>::quadraticNorm)
      .def("reverse", &Mvec<double>::reverse)
      .def("display", &Mvec<double>::display,
        py::call_guard<py::scoped_ostream_redirect,
                       py::scoped_estream_redirect>());

project_singular_metric_comment_begin
  mvec.def("outer_primal_dual", &Mvec<double>::outerPrimalDual);
  mvec.def("outer_dual_primal", &Mvec<double>::outerDualPrimal);
  mvec.def("outer_dual_dual", &Mvec<double>::outerDualDual);
  mvec.def("dual", &Mvec<double>::dual);
project_singular_metric_comment_end

  mvec.def("scalar_product", &Mvec<double>::scalarProduct);
  mvec.def("hestenes_product", &Mvec<double>::hestenesProduct);
  mvec.def("inv", &Mvec<double>::inv);

  mvec.def("grades", &Mvec<double>::grades);
  mvec.def("grade", [](const Mvec<double>& a){return a.grade();});
  mvec.def("grade", [](const Mvec<double>& a, const int i){return a.grade(i);});
  mvec.def("clear", &Mvec<double>::clear);

}

}  // namespace project_namespace
