// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// MetricTools.hpp
// This file is part of the Garamon Generator.
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

/// \file MetricTools.hpp
/// \author Stephane Breuils, Vincent Nozick
/// \brief Analyse and transform the metric related to the specified algebra.



#ifndef GARAGEN_METRICTOOLS_HPP
#define GARAGEN_METRICTOOLS_HPP

#include <Eigen/Core>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/LU> // for inverting

#include "Utilities.hpp"

bool isMatrixDiagonal(const Eigen::MatrixXd &metric, const double epsilon);

bool isMatrixIdentity(const Eigen::MatrixXd &metric, const double epsilon);

bool isMatrixPermutationOfDiagonal(const Eigen::MatrixXd &metric, const double epsilon);

unsigned int getRank(const Eigen::MatrixXd &metric);

void eigenDecomposition(const Eigen::MatrixXd &M, Eigen::MatrixXd &P, Eigen::MatrixXd &A);

double minAbsNonZeroValue(Eigen::VectorXd x);

Eigen::MatrixXd eigenRefinement(Eigen::MatrixXd &P, Eigen::MatrixXd &D, Eigen::MatrixXd &Pinv);

Eigen::MatrixXd numericalCleanUp(const Eigen::MatrixXd &M, const double epsilon);

Eigen::VectorXd vectorNumericalCleanUp(const Eigen::VectorXd& original, const double epsilon);

bool checkNumericalCleanUp(const Eigen::MatrixXd &M, const Eigen::MatrixXd &P, const Eigen::MatrixXd &A, const double epsilon);

bool checkNumericalCleanUp(const Eigen::MatrixXd &M, const Eigen::MatrixXd &P, const Eigen::MatrixXd &A, const Eigen::MatrixXd &Pinv, const double epsilon);

Eigen::SparseMatrix<double> numericalCleanUpSparse(const Eigen::MatrixXd &M, const double epsilon);

Eigen::SparseMatrix<double, Eigen::ColMajor> computePerGradeTransformationMatrix(const Eigen::MatrixXd &vectorTransformationMatrix,
                                         const unsigned int dimension, const unsigned int grade, const double epsilon);

Eigen::SparseMatrix<double, Eigen::ColMajor>  computeInverseTransformationMatrix(const Eigen::SparseMatrix<double, Eigen::ColMajor>& transformationMatrix, const double epsilon);

std::pair<std::vector<double>,std::vector<double>> computeTransformationMatricesToVector(const Eigen::MatrixXd &P, const double epsilon,std::vector<unsigned int>& transformationMatricesSizes,
                                                                                       std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor> >& allTransformationMatrices,
                                                                                       std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor> >& allInverseTransformationMatrices);

std::vector<double> transformationMatricesToVectorOfComponents(const Eigen::SparseMatrix<double, Eigen::ColMajor> &transformationMatrix, const int grade, const bool isInverse);

#endif //GARAGEN_METRICTOOLS_HPP
