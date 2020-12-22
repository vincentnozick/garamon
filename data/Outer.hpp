// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// Outer.hpp
// This file is part of the Garamon for project_namespace.
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

/// \file Outer.hpp
/// \author Stephane Breuils, Vincent Nozick
/// \brief Recursive outer product.


#ifndef project_inclusion_guard
#define project_inclusion_guard
#pragma once


#include "project_namespace/Mvec.hpp"

/*!
 * @namespace project_namespace
 */
namespace project_namespace {
    template<typename T> class Mvec;

    /// \brief Recursively compute the outer product between two multivectors mv1 and mv2, the result is put into the multivector mv3
    /// \tparam the type of value that we manipulate, either float or double or something.
    /// \param mv1 - the first multivector
    /// \param mv2 - the second multivector
    /// \param mv3 - the multivector that will content the result of the operation mv3 = mv1 ^ mv2
    /// \param grade_mv1 - the grade of the first multivector
    /// \param grade_mv2 - the grade of the second multivector
    /// \param grade_mv3 - the grade of the result which corresponds to grade_mv2+grade_mv1
    /// \param currentGradeMv1 - the current grade of the traversed tree of mv1
    /// \param currentGradeMv2 - the current grade of the traversed tree of mv2
    /// \param currentGradeMv3 - the current grade of the traversed tree of mv3
    /// \param sign - compute the sign of the outer product between two blades
    /// \param complement - activate the flip of sign
    /// \param indexLastVector_mv1 - last vector traversed in the multivector mv1
    /// \param indexLastVector_mv2 - last vector traversed in the multivector mv2
    /// \param indexLastVector_mv3 - last vector traversed in the multivector mv3
    template<typename T>
    void outerProductHomogeneous(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3, // homogeneous multivectors to be processed
                                 const unsigned int grade_mv1,const unsigned int grade_mv2,const unsigned int grade_mv3,         // grade of the k-vectors
                                 unsigned int currentXorIdx1=0, unsigned int currentXorIdx2=0, unsigned int currentXorIdx3=0,    // position in the prefix tree
                                 unsigned int currentGradeMv1=0, unsigned int currentGradeMv2=0, unsigned int currentGradeMv3=0, // grade relative to the position in the prefix tree
                                 int sign=1, int complement=1,
                                 unsigned int indexLastVector_mv1=0, unsigned int indexLastVector_mv2=0, unsigned int indexLastVector_mv3=1) {
        // sign updating
        int tmpSign = sign;
        if(complement == -1) {
            tmpSign = -tmpSign;
        }

        // if position in the tree for mv3 is not yet of grade of mv3, just call the recursive calls, without computation
        if(currentGradeMv3 < grade_mv3 ){
            // for each child of the node mv3 (call the recursive function) from the current node until the node has at least one child whose grade is grade_mv3
            for(unsigned int i = indexLastVector_mv3; (i<<(grade_mv3-(currentGradeMv3+1)))<(1<<algebraDimension);i*=2){

                unsigned int xorIndexMv1Child = currentXorIdx1 + i; // xor index of the child of the first tree multivector
                unsigned int xorIndexMv2Child = currentXorIdx2 + i; // xor index of the child of the second tree multivector
                unsigned int xorIndexMv3Child = currentXorIdx3 + i; // xor index of the child of the third tree multivector

                // if we do not reach the grade of mv1 AND if the child of the node of mv1 lead to at least one node whose grade is grade_mv1
                if ((currentGradeMv1 < grade_mv1) && (xorIndexMv1Child < (1 << algebraDimension))) {
                    outerProductHomogeneous<T>(mv1, mv2, mv3,
                                               grade_mv1, grade_mv2, grade_mv3,
                                               xorIndexMv1Child, currentXorIdx2, xorIndexMv3Child,
                                               currentGradeMv1 + 1, currentGradeMv2, currentGradeMv3 + 1,
                                               tmpSign, complement,
                                               i, indexLastVector_mv2, i << 1); //indexLastNodeTraversed_mv1+1,indexLastNodeTraversed_mv2
                }
                // if we do not reach the grade of mv2 AND if the child of the node of mv2 lead to at least one node whose grade is grade_mv2
                if ((currentGradeMv2 < grade_mv2) && (xorIndexMv2Child < (1 << algebraDimension))) {
                    outerProductHomogeneous<T>(mv1, mv2, mv3,
                                               grade_mv1, grade_mv2, grade_mv3,
                                               currentXorIdx1, xorIndexMv2Child, xorIndexMv3Child,
                                               currentGradeMv1, currentGradeMv2 + 1, currentGradeMv3 + 1,
                                               sign, -complement,
                                               indexLastVector_mv1, i, i << 1); //indexLastNodeTraversed_mv1, indexLastNodeTraversed_mv2+1
                }
            }
        } else { // currentGradeMv3 = grade_mv3 : compute the result and do not call the recursive process
            mv3(xorIndexToHomogeneousIndex[currentXorIdx3]) += sign * mv1(xorIndexToHomogeneousIndex[currentXorIdx1]) * mv2(xorIndexToHomogeneousIndex[currentXorIdx2]);
        }
    }

project_singular_metric_comment_begin
    /// \brief Recursively compute the outer product between two multivectors mv1 and mv2, we consider the dual form of the second multivector, the result is put into the multivector mv3
    /// \tparam the type of value that we manipulate, either float or double or something.
    /// \param mv1 - the first multivector
    /// \param mv2 - the second multivector as its dual form
    /// \param mv3 - the multivector that will contain the result of the operation mv3 = mv1 ^ dual(mv2)
    /// \param grade_mv1 - the grade of the first multivector
    /// \param grade_mv2 - the grade of the second multivector
    /// \param grade_mv3 - the grade of the result which corresponds to (dimension - grade_mv2) + grade_mv1
    /// \param currentGradeMv1 - the current grade of the traversed tree of mv1
    /// \param currentGradeMv2 - the current grade of the traversed tree of mv2
    /// \param currentGradeMv3 - the current grade of the traversed tree of mv3
    /// \param sign - compute the sign of the outer product between two blades
    /// \param complement - activate the flip of sign
    /// \param indexLastVector_mv1 - last vector traversed in the multivector mv1
    /// \param indexLastVector_mv2 - last vector traversed in the multivector mv2
    /// \param indexLastVector_mv3 - last vector traversed in the multivector mv3
    template<typename T>
    void outerProductPrimalDual(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3, // homogeneous multivectors to be processed
                                const unsigned int grade_mv1, const unsigned int grade_mv2, const unsigned int grade_mv3,   // grade of the k-vectors
                                unsigned int currentXorIdx1=0, unsigned int currentXorIdx2=0, unsigned int currentXorIdx3=(1<<algebraDimension)-1,  // position in the prefix tree
                                unsigned int currentGradeMv1=0, unsigned int currentGradeMv2=0, unsigned int currentGradeMv3=0, // grade relative to the position in the prefix tree
                                int sign=1, int complement=1,
                                unsigned int indexLastVector_mv1=0, unsigned int indexLastVector_mv2=1, unsigned int indexLastVector_mv3=0,
                                unsigned int depth=0) {
        // sign updating
        int tmpSign = sign;
        if(complement == -1) {
            tmpSign = -tmpSign;
        }

        // if position in the tree for mv3 is not yet of grade of mv3, just call the recursive calls, without computation
        if(currentGradeMv2 < grade_mv2 ){
            // for each child of the node mv3 (call the recursive function) from the current node until the node has at least one child whose grade is grade_mv3
            for(unsigned int i = indexLastVector_mv2; (i<<(grade_mv2-(currentGradeMv2+1)))<(1<<algebraDimension);i*=2){


                unsigned int xorIndexMv1Child = currentXorIdx1 + i;
                // the child of mv2 is not computed as in the wedge between two primal forms
                unsigned int xorIndexMv2Child = currentXorIdx2 + i; // e.g. e234 whose xor index is 1110b has e34 as first child which corresponds to 1100b = 1110b - 0010b
                // we traverse the multivector mv3 as before so the index of the child remains the same
                unsigned int xorIndexMv3Child = currentXorIdx3 - i;

                // if we do not reach the grade of mv1 AND if the child of the node of mv1 lead to at least one node whose grade is grade_mv1
                if ((currentGradeMv1 < grade_mv1) &&  ((i << (grade_mv1 - (currentGradeMv1 + 1))) < (1 << algebraDimension))) {
                    outerProductPrimalDual<T>(mv1, mv2, mv3,grade_mv1, grade_mv2, grade_mv3,
                                              xorIndexMv1Child, xorIndexMv2Child, currentXorIdx3,
                                              currentGradeMv1 + 1, currentGradeMv2 + 1, currentGradeMv3,
                                              tmpSign, complement,
                                              i, i << 1, indexLastVector_mv2,
                                              depth + 1); //indexLastNodeTraversed_mv1+1,indexLastNodeTraversed_mv2
                }
                // if we do not reach the grade of mv2 AND if the child of the node of mv2 lead to at least one node whose grade is grade_mv2
                if ((currentGradeMv3 < grade_mv3) && ((i << (grade_mv3 - (currentGradeMv3 + 1))) < (1 << algebraDimension)) ) {
                    outerProductPrimalDual<T>(mv1, mv2, mv3,grade_mv1, grade_mv2, grade_mv3,
                                              currentXorIdx1, xorIndexMv2Child, xorIndexMv3Child,
                                              currentGradeMv1, currentGradeMv2 + 1, currentGradeMv3 + 1,
                                              sign, -complement,
                                              indexLastVector_mv1, i << 1, i,
                                              depth+1); //indexLastNodeTraversed_mv1, indexLastNodeTraversed_mv2+1
                }
                tmpSign = -tmpSign;
                depth++;
            }
        } else { // currentGradeMv3 = grade_mv3 : compute the result and do not call the recursive process
            mv3(xorIndexToHomogeneousIndex[currentXorIdx3]) += recursiveDualCoefficients<T>[currentXorIdx2] * sign * mv1(xorIndexToHomogeneousIndex[currentXorIdx1]) * mv2(xorIndexToHomogeneousIndex[currentXorIdx2]);
        }
    }
    


    /// \brief Recursively compute the outer product between two multivectors mv1 and mv2, we consider the dual form of the first multivector, the result is put into the multivector mv3
    /// \tparam the type of value that we manipulate, either float or double or something.
    /// \param mv1 - the first multivector as its dual form
    /// \param mv2 - the second multivector 
    /// \param mv3 - the multivector that will contain the result of the operation mv3 = dual(mv1) ^ mv2
    /// \param grade_mv1 - the grade of the first multivector
    /// \param grade_mv2 - the grade of the second multivector
    /// \param grade_mv3 - the grade of the result which corresponds to (dimension - grade_mv2) + grade_mv1
    /// \param currentGradeMv1 - the current grade of the traversed tree of mv1
    /// \param currentGradeMv2 - the current grade of the traversed tree of mv2
    /// \param currentGradeMv3 - the current grade of the traversed tree of mv3
    /// \param sign - compute the sign of the outer product between two blades
    /// \param complement - activate the flip of sign
    /// \param indexLastVector_mv1 - last vector traversed in the multivector mv1
    /// \param indexLastVector_mv2 - last vector traversed in the multivector mv2
    /// \param indexLastVector_mv3 - last vector traversed in the multivector mv3
    template<typename T>
    void outerProductDualPrimal(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3, // homogeneous multivectors to be processed
                                const unsigned int grade_mv1, const unsigned int grade_mv2, const unsigned int grade_mv3,   // grade of the k-vectors
                                unsigned int currentXorIdx1=0, unsigned int currentXorIdx2=0, unsigned int currentXorIdx3=(1<<algebraDimension)-1,             // position in the prefix tree
                                unsigned int currentGradeMv1=0, unsigned int currentGradeMv2=0, unsigned int currentGradeMv3=0, // grade relative to the position in the prefix tree
                                int sign=1, int complement=1,
                                unsigned int indexLastVector_mv1=0, unsigned int indexLastVector_mv2=1, unsigned int indexLastVector_mv3=0,
                                unsigned int depth=0) {
        // sign updating
        int tmpSign = sign;
        if(complement == -1) {
            tmpSign = -tmpSign;
        }

        // if position in the tree for mv3 is not yet of grade of mv3, just call the recursive calls, without computation
        if(currentGradeMv2 < grade_mv2 ){
            // for each child of the node mv3 (call the recursive function) from the current node until the node has at least one child whose grade is grade_mv3
            for(unsigned int i = indexLastVector_mv2; (i<<(grade_mv2-(currentGradeMv2+1)))<(1<<algebraDimension);i*=2){

                unsigned int xorIndexMv1Child = currentXorIdx1 + i;
                // the child of mv2 is not computed as in the wedge between two primal forms
                unsigned int xorIndexMv2Child = currentXorIdx2 + i; // e.g. e234 whose xor index is 1110b has e34 as first child which corresponds to 1100b = 1110b - 0010b
                // we traverse the multivector mv3 as before so the index of the child remains the same
                unsigned int xorIndexMv3Child = currentXorIdx3 - i;

                // if we do not reach the grade of mv1 AND if the child of the node of mv1 lead to at least one node whose grade is grade_mv1
                if ((currentGradeMv1 < grade_mv1) &&  ((i << (grade_mv1 - (currentGradeMv1 + 1))) < (1 << algebraDimension))) {
                    outerProductDualPrimal<T>(mv1, mv2, mv3,grade_mv1, grade_mv2, grade_mv3,
                                              xorIndexMv1Child, xorIndexMv2Child, currentXorIdx3,
                                              currentGradeMv1 + 1, currentGradeMv2 + 1, currentGradeMv3,
                                              tmpSign, complement,
                                              i, i << 1, indexLastVector_mv2,
                                              depth + 1); //indexLastNodeTraversed_mv1+1,indexLastNodeTraversed_mv2
                }
                // if we do not reach the grade of mv2 AND if the child of the node of mv2 lead to at least one node whose grade is grade_mv2
                if ((currentGradeMv3 < grade_mv3) && ((i << (grade_mv3 - (currentGradeMv3 + 1))) < (1 << algebraDimension)) ) {
                    outerProductDualPrimal<T>(mv1, mv2, mv3,grade_mv1, grade_mv2, grade_mv3,
                                              currentXorIdx1, xorIndexMv2Child, xorIndexMv3Child,
                                              currentGradeMv1, currentGradeMv2 + 1, currentGradeMv3 + 1,
                                              sign, -complement,
                                              indexLastVector_mv1, i << 1, i,
                                              depth+1); //indexLastNodeTraversed_mv1, indexLastNodeTraversed_mv2+1
                }
                tmpSign = -tmpSign;
                depth++;
            }
        } else { // currentGradeMv3 = grade_mv3 : compute the result and do not call the recursive process
            mv3(xorIndexToHomogeneousIndex[currentXorIdx3]) += recursiveDualCoefficients<T>[currentXorIdx2] * sign * mv1(xorIndexToHomogeneousIndex[currentXorIdx1]) * mv2(xorIndexToHomogeneousIndex[currentXorIdx2]);
        }
    }

    /// \brief Recursively compute the outer product between two multivectors mv1 and mv2, we consider the dual form of the two multivectors, the result is put into the multivector mv3
    /// \tparam the type of value that we manipulate, either float or double or something.
    /// \param mv1 - the first multivector as its dual form
    /// \param mv2 - the second multivector as its dual form
    /// \param mv3 - the multivector that will contain the result of the operation mv3 = dual(mv1) ^ dual(mv2)
    /// \param grade_mv1 - the grade of the first multivector
    /// \param grade_mv2 - the grade of the second multivector
    /// \param grade_mv3 - the grade of the result which corresponds to (dimension - grade_mv2) + grade_mv1
    /// \param currentGradeMv1 - the current grade of the traversed tree of mv1
    /// \param currentGradeMv2 - the current grade of the traversed tree of mv2
    /// \param currentGradeMv3 - the current grade of the traversed tree of mv3
    /// \param sign - compute the sign of the outer product between two blades
    /// \param complement - activate the flip of sign
    /// \param indexLastVector_mv1 - last vector traversed in the multivector mv1
    /// \param indexLastVector_mv2 - last vector traversed in the multivector mv2
    /// \param indexLastVector_mv3 - last vector traversed in the multivector mv3
    template<typename T>
    void outerProductDualDual(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3, // homogeneous multivectors to be processed
                                const unsigned int grade_mv1, const unsigned int grade_mv2, const unsigned int grade_mv3,   // grade of the k-vectors
                                unsigned int currentXorIdx1=(1<<algebraDimension)-1, unsigned int currentXorIdx2=(1<<algebraDimension)-1, unsigned int currentXorIdx3=(1<<algebraDimension)-1,             // position in the prefix tree
                                unsigned int currentGradeMv1=0, unsigned int currentGradeMv2=0, unsigned int currentGradeMv3=0, // grade relative to the position in the prefix tree
                                int sign=1, int complement=1,
                                unsigned int indexLastVector_mv1=0, unsigned int indexLastVector_mv2=1, unsigned int indexLastVector_mv3=0,
                                unsigned int depth=0) {
        // sign updating
        int tmpSign = sign;
        if(complement == -1) {
            tmpSign = -tmpSign;
        }

        // if position in the tree for mv3 is not yet of grade of mv3, just call the recursive calls, without computation
        if(currentGradeMv2 < grade_mv2 ){
            // for each child of the node mv3 (call the recursive function) from the current node until the node has at least one child whose grade is grade_mv3
            for(unsigned int i = indexLastVector_mv2; (i<<(grade_mv2-(currentGradeMv2+1)))<(1<<algebraDimension);i*=2){

                unsigned int xorIndexMv1Child = currentXorIdx1 - i;
                // the child of mv2 is not computed as in the wedge between two primal forms
                unsigned int xorIndexMv2Child = currentXorIdx2 - i; // e.g. e234 whose xor index is 1110b has e34 as first child which corresponds to 1100b = 1110b - 0010b
                // we traverse the multivector mv3 as before so the index of the child remains the same
                unsigned int xorIndexMv3Child = currentXorIdx3 - i;

                // if we do not reach the grade of mv1 AND if the child of the node of mv1 lead to at least one node whose grade is grade_mv1
                if ((currentGradeMv1 < grade_mv1) &&  ((i << (grade_mv1 - (currentGradeMv1 + 1))) < (1 << algebraDimension))) {
                    outerProductDualDual<T>(mv1, mv2, mv3,grade_mv1, grade_mv2, grade_mv3,
                                              xorIndexMv1Child, xorIndexMv2Child, currentXorIdx3,
                                              currentGradeMv1 + 1, currentGradeMv2 + 1, currentGradeMv3,
                                              tmpSign, complement,
                                              i, i << 1, indexLastVector_mv2,
                                              depth + 1); //indexLastNodeTraversed_mv1+1,indexLastNodeTraversed_mv2
                }
                // if we do not reach the grade of mv2 AND if the child of the node of mv2 lead to at least one node whose grade is grade_mv2
                if ((currentGradeMv3 < grade_mv3) && ((i << (grade_mv3 - (currentGradeMv3 + 1))) < (1 << algebraDimension)) ) {
                    outerProductDualDual<T>(mv1, mv2, mv3,grade_mv1, grade_mv2, grade_mv3,
                                              currentXorIdx1, xorIndexMv2Child, xorIndexMv3Child,
                                              currentGradeMv1, currentGradeMv2 + 1, currentGradeMv3 + 1,
                                              sign, -complement,
                                              indexLastVector_mv1, i << 1, i,
                                              depth+1); //indexLastNodeTraversed_mv1, indexLastNodeTraversed_mv2+1
                }
                tmpSign = -tmpSign;
                depth++;
            }
        } else { // currentGradeMv3 = grade_mv3 : compute the result and do not call the recursive process
            mv3(xorIndexToHomogeneousIndex[currentXorIdx3]) += recursiveDualCoefficients<T>[currentXorIdx2] * sign * mv1(xorIndexToHomogeneousIndex[currentXorIdx1]) * mv2(xorIndexToHomogeneousIndex[currentXorIdx2]);
        }
    }


project_singular_metric_comment_end

}/// End of Namespace

#endif // project_inclusion_guard