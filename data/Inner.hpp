// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// Inner.hpp
// This file is part of the Garamon for project_namespace.
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

/// \file Inner.hpp
/// \author Stephane Breuils, Vincent Nozick
/// \brief Recursive inner product for project_namespace.

#ifndef project_inclusion_guard
#define project_inclusion_guard
#pragma once

#include "project_namespace/Mvec.hpp"
#include "project_namespace/Constants.hpp"

/*!
 * @namespace project_namespace
 */
namespace project_namespace {

    template<typename T> class Mvec;


    /// \brief Recursively compute the left contraction product between two multivectors mv1 and mv2, the result is put into the multivector mv3
    /// \tparam the type of value that we manipulate, either float or double or something.
    /// \param mv1 - the first multivector
    /// \param mv2 - the second multivector
    /// \param mv3 - the multivector that will content the result of the operation mv3 = mv1 * mv2
    /// \param grade_mv1 - the grade of the first multivector
    /// \param grade_mv2 - the grade of the second multivector
    /// \param grade_mv3 - the grade of the result which corresponds to grade_mv2-grade_mv1
    /// \param currentGradeMv1 - the current grade of the traversed tree of mv1
    /// \param currentGradeMv2 - the current grade of the traversed tree of mv2
    /// \param currentGradeMv3 - the current grade of the traversed tree of mv3
    /// \param sign - compute the sign of the geometric product between two blades
    /// \param complement - activate the flip of sign
    /// \param indexLastVector_mv1 - last vector traversed in the multivector mv1
    /// \param indexLastVector_mv2 - last vector traversed in the multivector mv2
    /// \param indexLastVector_mv3 - last vector traversed in the multivector mv3
    /// \param currentMetricCoefficient - coefficient that is used for handling the metric in the recursive formula
    /// \param depth - depth in the resulting multivector tree
    template<typename T>
    void leftContractionProductHomogeneous( const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3, // homogeneous multivectors to be processed
                                            const unsigned int grade_mv1, const unsigned int grade_mv2, const unsigned int grade_mv3,   // grade of the k-vectors
                                            unsigned int currentXorIdx1=0, unsigned int currentXorIdx2=0, unsigned int currentXorIdx3=0,             // position in the prefix tree
                                            unsigned int currentGradeMv1=0, unsigned int currentGradeMv2=0, unsigned int currentGradeMv3=0, // grade relative to the position in the prefix tree
                                            int sign=1, int complement=1,
                                            unsigned int indexLastVector_mv1=0, unsigned int indexLastVector_mv2=1, unsigned int indexLastVector_mv3=0, double currentMetricCoefficient=1.0, int depth=0) {
        // sign updating
        int tmpSign = sign;
        if(complement == -1) {
            tmpSign = -tmpSign;
        }

        // if position in the tree for mv3 is not yet of grade of mv2, just call the recursive calls, without computation
        if(currentGradeMv2 < grade_mv2 ){
            // for each child of the node mv2 (call the recursive function) from the current node until the node has at least one child whose grade is grade_mv2
            for(unsigned int i = indexLastVector_mv2; (i<<(grade_mv2-(currentGradeMv2+1)))<(1<<algebraDimension);i*=2){

                unsigned int xorIndexMv1Child = currentXorIdx1 + i; // xor index of the child of the first tree multivector
                unsigned int xorIndexMv2Child = currentXorIdx2 + i; // xor index of the child of the second tree multivector
                unsigned int xorIndexMv3Child = currentXorIdx3 + i; // xor index of the child of the third tree multivector

                // if we do not reach the grade of mv1 AND if the child of the node of mv1 lead to at least one node whose grade is grade_mv1
                if ( (currentGradeMv1 < grade_mv1) &&  ((i << (grade_mv1 - (currentGradeMv1 + 1))) < (1 << algebraDimension)) ) {
                    leftContractionProductHomogeneous<T>(mv1, mv2, mv3,
                                               grade_mv1, grade_mv2, grade_mv3,
                                               xorIndexMv1Child, xorIndexMv2Child, currentXorIdx3,
                                               currentGradeMv1 + 1, currentGradeMv2 +1 , currentGradeMv3 ,
                                               tmpSign, -complement,
                                               i, i << 1, indexLastVector_mv2,
                                                         diagonalMetric<T>(depth)*currentMetricCoefficient, depth+1);
                }

                // if we do not reach the grade of mv3 AND if the child of the node of mv3 lead to at least one node whose grade is grade_mv3
                if ((currentGradeMv3 < grade_mv3) && ((i << (grade_mv3 - (currentGradeMv3 + 1))) < (1 << algebraDimension)) ) {
                    leftContractionProductHomogeneous<T>(mv1, mv2, mv3,
                                               grade_mv1, grade_mv2, grade_mv3,
                                               currentXorIdx1, xorIndexMv2Child, xorIndexMv3Child,
                                               currentGradeMv1, currentGradeMv2 + 1, currentGradeMv3 + 1,
                                               sign, -complement,
                                               indexLastVector_mv1, i << 1, i,
                                               currentMetricCoefficient, depth+1);
                }
                depth++;
            }
        } else { // currentGradeMv2 = grade_mv2 : compute the result and do not call the recursive process
            mv3(xorIndexToHomogeneousIndex[currentXorIdx3]) += sign * currentMetricCoefficient * mv1(xorIndexToHomogeneousIndex[currentXorIdx1]) * mv2(xorIndexToHomogeneousIndex[currentXorIdx2]);
        }
    }



    /// \brief Recursively compute the right contraction product between two multivectors mv1 and mv2, the result is put into the multivector mv3
    /// \tparam the type of value that we manipulate, either float or double or something.
    /// \param mv1 - the first multivector
    /// \param mv2 - the second multivector
    /// \param mv3 - the multivector that will content the result of the operation mv3 = mv1 * mv2
    /// \param grade_mv1 - the grade of the first multivector
    /// \param grade_mv2 - the grade of the second multivector
    /// \param grade_mv3 - the grade of the result which corresponds to grade_mv1-grade_mv2
    /// \param currentGradeMv1 - the current grade of the traversed tree of mv1
    /// \param currentGradeMv2 - the current grade of the traversed tree of mv2
    /// \param currentGradeMv3 - the current grade of the traversed tree of mv3
    /// \param sign - compute the sign of the geometric product between two blades
    /// \param complement - activate the flip of sign
    /// \param indexLastVector_mv1 - last vector traversed in the multivector mv1
    /// \param indexLastVector_mv2 - last vector traversed in the multivector mv2
    /// \param indexLastVector_mv3 - last vector traversed in the multivector mv3
    /// \param currentMetricCoefficient - coefficient that is used for handling the metric in the recursive formula
    /// \param depth - depth in the resulting multivector tree
    template<typename T>
    void rightContractionProductHomogeneous(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3, // homogeneous multivectors to be processed
                                           const unsigned int grade_mv1,const unsigned int grade_mv2,const unsigned int grade_mv3,   // grade of the k-vectors
                                           unsigned int currentXorIdx1=0, unsigned int currentXorIdx2=0, unsigned int currentXorIdx3=0,             // position in the prefix tree
                                           unsigned int currentGradeMv1=0, unsigned int currentGradeMv2=0, unsigned int currentGradeMv3=0, // grade relative to the position in the prefix tree
                                           int sign=1, int complement=1,
                                           unsigned int indexLastVector_mv1=1, unsigned int indexLastVector_mv2=0, unsigned int indexLastVector_mv3=0, double currentMetricCoefficient=1.0, int depth=0) {

        // sign updating
        int tmpSign = sign;
        if(complement == -1) {
            tmpSign = -tmpSign;
        }

        // if position in the tree for mv3 is not yet of grade of mv3, just call the recursive calls, without computation
        if(currentGradeMv1 < grade_mv1 ){
            // for each child of the node mv1 (call the recursive function) from the current node until the node has at least one child whose grade is grade_mv1
            for(unsigned int i = indexLastVector_mv1; (i<<(grade_mv1-(currentGradeMv1+1)))<(1<<algebraDimension);i*=2){
                unsigned int xorIndexMv1Child = currentXorIdx1 + i; // xor index of the child of the first tree multivector
                unsigned int xorIndexMv2Child = currentXorIdx2 + i; // xor index of the child of the second tree multivector
                unsigned int xorIndexMv3Child = currentXorIdx3 + i; // xor index of the child of the third tree multivector

                // if we do not reach the grade of mv2 AND if the child of the node of mv2 lead to at least one node whose grade is grade_mv2
                if ( (currentGradeMv2 < grade_mv2) &&  ((i << (grade_mv2 - (currentGradeMv2 + 1))) < (1 << algebraDimension)) ) {
                    rightContractionProductHomogeneous<T>(mv1, mv2, mv3,
                                                         grade_mv1, grade_mv2, grade_mv3,
                                                         xorIndexMv1Child, xorIndexMv2Child, currentXorIdx3,
                                                         currentGradeMv1 + 1, currentGradeMv2 +1 , currentGradeMv3 ,
                                                         tmpSign, -complement,
                                                          i << 1,i, indexLastVector_mv3,
                                                         diagonalMetric<T>(depth)*currentMetricCoefficient,depth+1);
                }

                // if we do not reach the grade of mv3 AND if the child of the node of mv3 lead to at least one node whose grade is grade_mv3
                if ((currentGradeMv3 < grade_mv3) && ((i << (grade_mv3 - (currentGradeMv3 + 1))) < (1 << algebraDimension)) ) {
                    rightContractionProductHomogeneous<T>(mv1, mv2, mv3,
                                                         grade_mv1, grade_mv2, grade_mv3,
                                                          xorIndexMv1Child, currentXorIdx2, xorIndexMv3Child,
                                                         currentGradeMv1+1, currentGradeMv2, currentGradeMv3 + 1,
                                                          tmpSign, complement,
                                                          i << 1, indexLastVector_mv2, i,
                                                          currentMetricCoefficient, depth+1);
                }
                depth++;
            }
        } else { // currentGradeMv1 = grade_mv1 : compute the result and do not call the recursive process
           mv3(xorIndexToHomogeneousIndex[currentXorIdx3]) += sign * currentMetricCoefficient * mv1(xorIndexToHomogeneousIndex[currentXorIdx1]) * mv2(xorIndexToHomogeneousIndex[currentXorIdx2]);
        }
    }



 
} /// End of Namespace

#endif // project_inclusion_guard
