// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// GeometricExplicit.hpp
// This file is part of the Garamon for project_namespace.
// Authors: Stephane Breuils and Vincent Nozick
// Conctat: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

/// \file GeometricExplicit.hpp
/// \author Stephane Breuils, Vincent Nozick
/// \brief Explicit precomputed per grades geometric products of project_namespace.


#ifndef project_inclusion_guard
#define project_inclusion_guard
#pragma once

#include <Eigen/Core>

#include "project_namespace/Mvec.hpp"
#include "project_namespace/Constants.hpp"


/*!
 * @namespace project_namespace
 */
namespace project_namespace {
    template<typename T> class Mvec;

    project_explicit_geometric_functions
    project_explicit_geometric_pointer_functions

}/// End of Namespace

#endif // project_inclusion_guard