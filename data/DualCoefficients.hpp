// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// DualCoefficients.hpp
// This file is part of the Garamon for project_namespace.
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

// \file DualCoefficients.hpp
// \brief this files generates and load the elements of the fast dual into array of Eigen matrices 
// \author V. Nozick, S. Breuils

#ifndef project_inclusion_guard
#define project_inclusion_guard
#pragma once


#include <Eigen/Core>
#include <vector>
#include <array>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <string>


/*!
 * @namespace project_namespace
 */
namespace project_namespace{

	/// Decode the string encodedPerGradeMatrixComponents to a vector of T. This vector will contain the coefficients used to compute the dual
	template<typename T>
	const std::vector<T> decodeStringOfDualCoefToVecOfT(std::string encodedPerGradeMatrixComponents) {
	    std::vector<T> resultDecodedGrade;
	    std::istringstream streamStringOverFloat(encodedPerGradeMatrixComponents);
	    // parse as double (full precision), then convert to T
	    std::copy(std::istream_iterator<double>(streamStringOverFloat),
		std::istream_iterator<double>(),
		std::back_inserter(resultDecodedGrade));
	    return resultDecodedGrade;
	}


	/// load the coefficients of the dual into array of Eigen matrices 
project_fill_dual_array

}/// End of Namespace

#endif // project_inclusion_guard
