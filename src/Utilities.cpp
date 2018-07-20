// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// Utilities.cpp
// This file is part of the Garamon Generator.
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

#include "Utilities.hpp"


// compute the factorial
unsigned int factorial(unsigned int n) {
    return (n <= 1) ? 1 : n * factorial(n - 1);
}


// New version: Iterative version of the binomial Coefficient without computing the factorial
unsigned int bin_coeff(unsigned int n, unsigned int k) {
    if(k>n) return 0;

    unsigned int res = 1;
    // Because C(n, k) = C(n, n-k)
    if (k > n - k)
        k = n - k;
    // Compute value of [n * (n-1) * ... * (n-k+1)] / [k * (k-1) *----* 1]
    for (unsigned int i = 0; i < k; ++i) {
        res *= (n - i);
        res /= (i + 1);
    }
    return res;
}


// use the binomial coefficient to comute
void computePerGradeStartingIndex(const unsigned int &dimension, std::vector<int> &perGradeStartingIndex,
                                  unsigned int currentIndex, unsigned int currentGrade){

    currentIndex += bin_coeff(dimension, currentGrade);
    perGradeStartingIndex.push_back(currentIndex);
    if(currentGrade < (dimension-1)){
        computePerGradeStartingIndex(dimension,perGradeStartingIndex,currentIndex, currentGrade+1);
    }
}


unsigned int hammingWeight(unsigned int word){
    return __builtin_popcount(word);
}


int computeSign(const unsigned int mv1, const unsigned int mv2){
    unsigned int nbOfOnes = 0;
    unsigned int tmpmva = mv1;
    tmpmva = tmpmva >> 1;
    while(tmpmva!=0){
        // Count number of ones in common
        nbOfOnes += hammingWeight(tmpmva & mv2);
        tmpmva = tmpmva >> 1;
    }
    return (-2*(nbOfOnes%2))+1;
}


// mv1 and mv2 are binary expression of some basis vectors (i.e. e23 -> 0110)
// this function computes the binary expression of the basis vector containing the result of mv1^mv2,
// as well as the sign corresponding to the number of permutations required to order mv1^mv2.
unsigned int outerProductBinaryComponents(const unsigned int mv1, const unsigned int mv2, int &sign){
    sign = computeSign(mv1, mv2);
    return mv1^mv2;
}


void getSignAndBladeOuterProduct(unsigned int& mvC, const unsigned int mvA, const unsigned int mvB, double& coeffC, const double coeffA, const double coeffB){

    int sign = 1;

    if( (mvA & mvB) == 0){ // if mvA and mvB has some common vector basis, the outer product will be 0
        mvC = outerProductBinaryComponents(mvA, mvB, sign);
        coeffC += sign * coeffA * coeffB;
    }
}


std::vector<unsigned int> getSetOfCombinationsFromXorIndexation(const int dim, const std::vector<std::vector<unsigned int> >& combinations){

    std::vector<unsigned int> setOfCombinationsFromXorIndex((unsigned int)1 << dim);

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
