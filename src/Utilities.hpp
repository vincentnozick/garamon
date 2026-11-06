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

#include <string>
#include <vector>


// Iterative version of the binomial Coefficient without computing the factorial
unsigned int bin_coeff(unsigned int n, unsigned int k);

// use the binomial coefficient to comute
void computePerGradeStartingIndex(const unsigned int &dimension, std::vector<int> &perGradeStartingIndex,
                                  unsigned int currentIndex, unsigned int currentGrade);

// (the blade products, e.g. the outer product sign, are in ProductTools.hpp)

// combinations is the sequence of all possible elements of a considered grade
// i.e. for grade 2, dim 3 :                (1,2) (1,3) (2,3)
// order induced by the xor:                scal, 1, 2, 12, 3, 13, 23, 123
// return their position in the sequence:      -, -, -,  0, -,  1,  2, -
std::vector<unsigned int> getSetOfCombinationsFromXorIndexation(const int dim, const std::vector<std::vector<unsigned int> >& combinations);

// generate the sequence of all possible elements of a considered grade
// e.g. generateCombinations(4,2) will result to the sequence (1,2), (1,3), (1,4), (2,3), (2,4), (3,4)
// here 'n' is the dimension and 'r' is the considered grade
std::vector<std::vector<unsigned int> > generateCombinations(const unsigned int n, const unsigned int r);



// convert a double to a string without any precision loss (unlike std::to_string that keeps only 6 decimals), i.e. "0.5", "-2", "0.33333333333333331"
std::string doubleToString(const double value);

#endif //GARAMON_UTILITY_HPP_H
