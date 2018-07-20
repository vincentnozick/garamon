// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// OuterExplicit.hpp
// This file is part of the Garamon for project_namespace.
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

/// \file OuterExplicit.hpp
/// \author Stephane Breuils, Vincent Nozick
/// \brief Explicit precomputed per grades outer product.


#ifndef project_inclusion_guard
#define project_inclusion_guard
#pragma once

#include <Eigen/Core>

#include "project_namespace/Mvec.hpp"
#include "project_namespace/Outer.hpp"


/*!
 * @namespace project_namespace
 */
namespace project_namespace {
    template<typename T> class Mvec;

    project_explicit_outer_functions

    project_explicit_outer_pointer_functions

}/// End of Namespace

#endif // project_inclusion_guard