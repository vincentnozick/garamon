// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// Utilities.cpp
// This file is part of the Garamon Generator.
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

#include "Utilities.hpp"

#include <algorithm>
#include <cstddef>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>


// New version: Iterative version of the binomial Coefficient without computing the factorial
unsigned int bin_coeff(unsigned int n, unsigned int k) {
    if(k>n) return 0;

    // 64 bits intermediate results: res * (n-i) overflows 32 bits from C(31,15) on
    unsigned long long res = 1;
    // Because C(n, k) = C(n, n-k)
    if (k > n - k)
        k = n - k;
    // Compute value of [n * (n-1) * ... * (n-k+1)] / [k * (k-1) *----* 1]
    for (unsigned int i = 0; i < k; ++i) {
        res *= (n - i);
        res /= (i + 1);
    }
    return (unsigned int)res;
}


// use the binomial coefficient to comute
void computePerGradeStartingIndex(const unsigned int &dimension, std::vector<int> &perGradeStartingIndex,
                                  unsigned int currentIndex, unsigned int currentGrade){

    currentIndex += bin_coeff(dimension, currentGrade);
    perGradeStartingIndex.push_back(currentIndex);
    if(currentGrade + 1 < dimension){ // (not 'currentGrade < dimension-1': unsigned underflow for dimension 0)
        computePerGradeStartingIndex(dimension,perGradeStartingIndex,currentIndex, currentGrade+1);
    }
}


std::vector<unsigned int> getSetOfCombinationsFromXorIndexation(const int dim, const std::vector<std::vector<unsigned int> >& combinations){

    std::vector<unsigned int> setOfCombinationsFromXorIndex(std::size_t(1) << dim);

    // for each possible sequence
    for(unsigned int i=0; i<combinations.size(); ++i){
        unsigned int idx=0;

        // for each element of a sequence
        for(unsigned int j=0; j<combinations[i].size(); ++j) {
            idx += 1 << (combinations[i][j]);
        }

        // the sequence defined by idx, the corresponding index in the sequence is i
        setOfCombinationsFromXorIndex[idx] = i;
    }

    return setOfCombinationsFromXorIndex;
}


std::vector<std::vector<unsigned int> > generateCombinations(const unsigned int n, const unsigned int r){

    // n is the dimension of the algebra
    std::vector<bool> v(n);
    std::fill(v.begin(), v.end(), false);

    // the r (the considered grade) first elements of the n elements of v are true
    // v represent the first element of the sequense (i.e (1,2) for r=2)
    std::fill(v.begin(), v.begin() + r, true);

    // output sequence of the each combinaison of r-vectors
    std::vector<std::vector<unsigned int> > sequence;

    do {
        // one element of the sequence is composed of the r elements of v
        std::vector<unsigned int> kvector;
        for(unsigned int i=0; i<n; ++i) {
            if(v[i]) {
                kvector.push_back((unsigned int)i);
            }
        }

        // add this element to the sequence
        sequence.push_back(kvector);

    } while(std::prev_permutation(v.begin(), v.end())); // compute next permutation of the true values of v

    return sequence;
}


std::string doubleToString(const double value){
    std::ostringstream stream;
    stream << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
    return stream.str();
}
