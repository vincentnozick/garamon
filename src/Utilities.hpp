// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// Utilities.hpp
// This file is part of the Garamon Generator.
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

/// \file Utilities.hpp
/// \author Stephane Breuils, Vincent Nozick
/// \brief Some basic and generic tools.


#ifndef GARAMON_UTILITY_HPP_H
#define GARAMON_UTILITY_HPP_H

#include "ProductTools.hpp"


// compute the factorial
unsigned int factorial(unsigned int n);

// Iterative version of the binomial Coefficient without computing the factorial
unsigned int bin_coeff(unsigned int n, unsigned int k);

// use the binomial coefficient to comute
void computePerGradeStartingIndex(const unsigned int &dimension, std::vector<int> &perGradeStartingIndex,
                                  unsigned int currentIndex, unsigned int currentGrade);

// TO REMOVE if redondant with a similiar function in ProductTools.cpp
// Compute the number of 1s in the binary integer word
unsigned int hammingWeight(unsigned int word);

// TO REMOVE if redondant with a similiar function in ProductTools.cpp
// mva = index of the basis vector used (equivalent to "each bit corresponds to a vector basis")
// for e23 => e2 and e3 will be set to 1, the rest to 0
// sign = (-1)^{grade(mva)*grade(mvb)/2}
int computeSign(const unsigned int mv1, const unsigned int mv2);

// TO REMOVE if redondant with a similiar function in ProductTools.cpp
// compute the sign and index of the result of the wedge product of mv1 and mv2
unsigned int outerProductBinaryComponents(const unsigned int mv1, const unsigned int mv2, int &sign);

// TO REMOVE if redondant with a similiar function in ProductTools.cpp
// same as perBladeOuterProduct, but also includes the coefficient of the result
void getSignAndBladeOuterProduct(unsigned int& mvC, const unsigned int mvA, const unsigned int mvB, double& coeffC, const double coeffA, const double coeffB);

// TO REMOVE if redondant with a similiar function in ProductTools.cpp
// combinations is the sequence of all possible elements of a considered grade
// i.e. for grade 2, dim 3 :                (1,2) (1,3) (2,3)
// order induced by the xor:                scal, 1, 2, 12, 3, 13, 23, 123
// return their position in the sequence:      -, -, -,  0, -,  1,  2, -
std::vector<unsigned int> getSetOfCombinationsFromXorIndexation(const int dim, const std::vector<std::vector<unsigned int> >& combinations);

// TO REMOVE if redondant with a similiar function in ProductTools.cpp
// generate the sequence of all possible elements of a considered grade
// e.g. generateCombinations(4,2) will result to the sequence (1,2), (1,3), (1,4), (2,3), (2,4), (3,4)
// here 'n' is the dimension and 'r' is the considered grade
std::vector<std::vector<unsigned int> > generateCombinations(const unsigned int n, const unsigned int r);



#endif //GARAMON_UTILITY_HPP_H
