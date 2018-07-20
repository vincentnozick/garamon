// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// Constants.hpp
// This file is part of the Garamon for project_namespace.
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

/// \file Constants.hpp
/// \author Stephane Breuils, Vincent Nozick
/// \brief Constant values and data related to the specified geometric algebra (project_namespace)


// Doxygen
/// \version 0.1
/// \mainpage
/// \tableofcontents
/// \section instroduction_sec What for?
/// Garamon is a C++ library to represent an manipulate the geometric algebra objects.
/// \section install_bigsec How to install
/// \subsection dependencies_sec Dependecies
/// \li cmake (at least version 3.10)
/// \li Eigen (at least version 3.3.4)
/// \li Doxygen (if you want the documentation)
/// \subsection install_sec Install with cmake (Linux / MacOs)
/// \li go to garamon dir
/// \li mkdir build
/// \li cd build
/// \li cmake ..
/// \li make
/// \li (optional) sudo make install
/// \li (optional for documentation) make html
/// \li The documentation is located in [path to build]/doc/doc-doxygen/html/index.html


#ifndef project_inclusion_guard
#define project_inclusion_guard
#pragma once


#include <array>
#include <vector>
#include <utility>
#include <string>
#include <Eigen/Sparse>

#include "project_namespace/Utility.hpp"
#include "project_namespace/BasisTransformations.hpp"
#include "project_namespace/DualCoefficients.hpp"


/*!
 * @namespace project_namespace
 */
namespace project_namespace{

    constexpr unsigned int algebraDimension = project_algebra_dimension; /*!< dimension of the algebra (number of  basis vectors of grade 1) */

    constexpr unsigned int project_per_grade_starting_index;  /*!< array specifying the index of each first element of grade k in the full multivector */

    constexpr unsigned int project_array_binomial_coefficient;  /*!< array of the (dimension + 1) first binomial coefficients */

project_array_xorIndexConversion

project_dual_arrays_permutations_and_coefficients
project_dual_arrays_recursive_coefficients

project_pseudo_scalar_inverse
    const int project_sign_reverse; /*!< array of signs to avoid the computation of (-1)^k*(k-1)/2 during the reverse operation */

    const std::vector<std::string> basisVectors = project_basis_vectors_string; /*!< name of the basis vectors (of grade 1) */

    const std::string metric =
project_metric; /*!< metric / quadratic form of the algebra (inner product between basis vectors) */

    const unsigned int scalar = 0;
project_basis_vector_index    /*!< defines the constants for the cga */

project_load_transformation_matrices

    project_diagonal_Metric


}  // namespace


#endif // project_inclusion_guard
