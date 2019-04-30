// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// ProductToString.hpp
// This file is part of the Garamon Generator.
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

/// \file ProductToString.hpp
/// \author Stephane Breuils, Vincent Nozick
/// \brief Tools to convert numerical data in c++ code.


#ifndef GARAGEN_PRODUCTTOSTRING_HPP
#define GARAGEN_PRODUCTTOSTRING_HPP

#include <string>
#include <vector>

#include "MetaData.hpp"
#include "ProductTools.hpp"


std::string basisVectorsToString(const MetaData &metaData);

std::string perGradeStartingIndexToString(const std::vector<int> &perGradeStartingIndex);

std::string vectorIdToString(const MetaData &metaData, const bool lowercase, const int vectorNumber);

std::string metricToString(const MetaData &metaData);

std::string binomialCoefToString(const unsigned int &dimension);


std::string initKVectorDirectTransformationMatrixFunctionPrototype(const unsigned int grade);

std::string initKVectorInverseTransformationMatrixFunctionPrototype(const unsigned int grade);

std::string perGradetransformMatricesToString(const std::vector<double>& transformComponents, const unsigned int grade,  const unsigned int startingPosition, const unsigned int endPosition);

std::string perGradeEigenSparseMatrixInitialization(const std::vector<double>& transformComponents, const unsigned int dimension, const unsigned int grade,  const unsigned int startingPosition, const unsigned int endPosition, const unsigned int numberOfNonNullComponents, const std::string& nameComponentsDataStructure);

std::string loadAllDirectOrInverseMatrices(const std::vector<unsigned int> &sizeTransformationMatrices,
                                           const std::vector<double> &transformComponents, const bool isInverse);

std::string callAllDirectOrInverseMatricesFunctions(const std::vector<unsigned int>& sizeTransformationMatrices, const bool isInverse);

std::string loadAllDualCoefficientsArray(const std::vector<int> & sizeDualCoefficientsArray, const std::string& stringDualCoefficientsComponents);

std::string dualPermutationToString(const std::vector<double>& basisChangesComponents);

std::string staticOneComponentMultivectorPrototypePython();

std::string staticOneComponentMultivectorPrototype();

std::string oneComponentMultivectorPrototype();

std::string constantsDefinition();

std::string singularMetricCommentBegin();

std::string singularMetricCommentEnd();

std::string reverseSignArrayToString(const unsigned int &dimension);

std::string xorIndexToGradeAndHomogeneousIndexArraysToString(unsigned int dimension, const ProductTools& product);

std::string fastDualUtilities(unsigned int dimension, const ProductTools& product,
                                        const Eigen::VectorXd &diagonalMetric, double scaleInversePseudoScalar,std::string srcDirectory,
                                        std::string& fastDualComponents);


std::string fastDualUtilitiesBasisChange(unsigned int dimension, const ProductTools& product,
                                         const std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor> >& transformationMatrices,
                                         const std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor> >& inverseTransformationMatrices,
                                         const Eigen::VectorXd &diagonalMetric,
                                         double scaleInversePseudoScalar,
                                         std::string srcDirectory,
                                         std::string& fastDualComponents);

std::string primalWedgeDualUtilities(const unsigned int dimension, const ProductTools& product,
                                     const Eigen::VectorXd &diagonalMetric,
                                     double scaleInversePseudoScalar);

std::string primalWedgeDualUtilitiesBasisChange(unsigned int dimension, const ProductTools& product,
                                     const std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor> >& transformationMatrices,
                                     const std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor> >& inverseTransformationMatrices,
                                     const Eigen::VectorXd &diagonalMetric,
                                     double scaleInversePseudoScalar);

std::string basisTransformMatricesLoad();

std::string multivectorComponentBuilder(const MetaData& metaData, const std::string &data);

std::string diagonalMetricToString(const MetaData& metaData);

std::string templateVectorToOrthogonalBasisToString();

std::string templateTransformVectorToOrthogonalBasisToString();

std::string templateTransformVectorToOriginalBasisToString();

std::string outerProductExplicitComments(const unsigned int gradeMv1, const unsigned int gradeMv2, const unsigned int gradeMv3);

std::string outerProductExplicitPrototype(const unsigned int gradeMv1, const unsigned int gradeMv2);

std::string productListToString(std::list<productComponent<double>> &listOfExplicitProducts);

std::string generateOuterRecursive(const unsigned int gradeMv1, const unsigned int gradeMv2, const unsigned int gradeMv3);

std::string generateOuterExplicit_cpp(const MetaData &metaData, const ProductTools& product);

std::string generateOuterExplicitFunctionsPointer(const unsigned int &dimension);

std::string generateOuterExplicit_simd_64(const MetaData &metaData, const ProductTools& product);

std::string generateOuterExplicitFunctionsPointer_simd_64(const unsigned int &dimension);

std::string generateOuterExplicit_simd_32(const MetaData &metaData, const ProductTools& product);

std::string generateOuterExplicitFunctionsPointer_simd_32(const unsigned int &dimension);

std::string innerProductExplicitComments(const unsigned int gradeMv1, const unsigned int gradeMv2, const unsigned int gradeMv3);

std::string generateInnerRecursiveBasisChange(const unsigned int gradeMv1, const unsigned int gradeMv2, const unsigned int gradeMv3,
                                              const unsigned int dimension);

std::string generateInnerRecursiveBasisChangeFloat64(const unsigned int gradeMv1, const unsigned int gradeMv2, const unsigned int gradeMv3,
                                              const unsigned int dimension);

std::string generateInnerRecursive(const unsigned int gradeMv1, const unsigned int gradeMv2, const unsigned int gradeMv3);

std::string generateInnerRecursiveFloat64(const unsigned int gradeMv1, const unsigned int gradeMv2, const unsigned int gradeMv3);

std::string generateInnerExplicit_cpp(const MetaData &metaData,
                                      const ProductTools& product);

std::string generateInnerExplicitBasisChange_cpp(const MetaData &metaData,
                                      std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor> > transformationMatrices,
                                      std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor> > inverseTransformationMatrices,
                                      const ProductTools& product);

std::string generateInnerExplicitFunctionsPointer(const unsigned int &dimension);

std::string generateInnerExplicit_simd_64(const MetaData &metaData,
                                          const ProductTools& product);

std::string generateInnerExplicitBasisChange_simd_64(const MetaData &metaData,
                                          std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor> > transformationMatrices,
                                          std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor> > inverseTransformationMatrices,
                                          const ProductTools& product);

std::string generateInnerExplicitFunctionsPointer_simd_64(const unsigned int &dimension);

std::string geometricProductExplicitComments(const unsigned int gradeMv1, const unsigned int gradeMv2, const unsigned int gradeMvResult);


std::string generateGeometricExplicitBasisChange_cpp(const MetaData &metaData,
                                          std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor>> transformationMatrices,
                                          std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor>> inverseTransformationMatrices,
                                          const ProductTools& product);

std::string generateGeometricExplicit_cpp(const MetaData &metaData,
                                          const ProductTools& product);

std::string generateGeometricExplicitFunctionsPointer(const MetaData &metaData);

std::string recursiveGeometricProductCallToString(const unsigned int maxSize, bool changeMetricToOrtrhogonal);

std::string recursiveGeometricProductCallToStringFloat64(const unsigned int maxSize, bool changeMetricToOrtrhogonal);

std::string generateGeometricExplicit_simd_64(const MetaData &metaData,
                                              const ProductTools& product);

std::string generateGeometricExplicitBasisChange_simd_64(const MetaData &metaData,
                                              std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor>> transformationMatrices,
                                              std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor>> inverseTransformationMatrices,
                                              const ProductTools& product);

std::string generateGeometricExplicitFunctionsPointer_simd_64(const MetaData &metaData);


#endif //GARAGEN_PRODUCTTOSTRING_HPP
