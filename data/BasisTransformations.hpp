// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// Inner.hpp
// This file is part of the Garamon for project_namespace.
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

// \file BasisTransformations.hpp
// \author V. Nozick, S. Breuils
// \brief this files generates and load the elements of the transformation matrices into array of Eigen matrices 


#ifndef project_inclusion_guard
#define project_inclusion_guard
#pragma once


#include <Eigen/Sparse>
#include <vector>
#include <array>
#include <algorithm>
#include <iterator>
#include <sstream>


/*!
 * @namespace project_namespace
 */
namespace project_namespace{

	/// Decode the string encodedPerGradeMatrixComponents to a vector of T. This vector will contain the components of the transformation matrices
	template<typename T>
	const std::vector<T> decodeStringToVecOfT(std::string encodedPerGradeMatrixComponents) {
	    std::vector<T> resultDecodedGrade;
	    std::istringstream streamStringOverFloat(encodedPerGradeMatrixComponents);
	    std::copy(std::istream_iterator<float>(streamStringOverFloat),
		std::istream_iterator<float>(),
		std::back_inserter(resultDecodedGrade));
	    return resultDecodedGrade;
	}

project_basischange_direct_loading

project_basischange_inverse_loading

	/// initialize all the direct transformation matrices using array of eigen sparse matrices
project_calltobasischange_direct_loading

	/// initialize all the inverse transformation matrices using array of eigen sparse matrices
project_calltobasischange_inverse_loading

}/// End of Namespace

#endif // project_inclusion_guard
