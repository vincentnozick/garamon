// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// MetaData.hpp
// This file is part of the Garamon Generator.
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

/// \file MetaData.hpp
/// \author Stephane Breuils, Vincent Nozick
/// \brief Contains all the required information to define and create an algebra (dimension, metric, name, code options, ...)


#ifndef GARAGEN_METADATA_HPP__
#define GARAGEN_METADATA_HPP__

#include <iostream>
#include <vector>
#include <string>
#include <Eigen/Core>
#include "MetricTools.hpp"

// data required to build a Geometric Algebra Library
class MetaData
{
public:

    // namespace of the generated lib (ex: cga::)
    std::string namespaceName;

    // dimension of the algebra vector space (grade 1)
    unsigned int dimension;

    // metric = inner product of vectors, should be a symetric matrix
    Eigen::MatrixXd metric;

    // name of each basis vector, without the first letter (1,2,3 => e1, e2, e3, e12, ...).
    // for hig dimensions, add a '_' at the end of the name to avoid: 1, 2, 3, ..., 16 => e114 : 11,4 ou 1,14 ?? => e11_4_
    std::vector<std::string> basisVectorName;

    // true if the proposed metric is diagonal
    bool inputMetricDiagonal;

    // true if the final (or initial) diagonal metric is actually identity
    bool identityMetric;

    // true if the metric has rank = dimension of the vector space supporting the algebra
    bool fullRankMetric;

    // true if the initial metric is a permutation of a diagonal matrix
    bool inputMetricPermutationOfDiagonal;

    // numerical refinement of the eigen vectors / values
    bool useEigenRefinement;

    // replace near zeros by zeros, near integers by integers, ...
    bool useNumericalCleanUp;

    // Each multivector component of grade k can have a dedicated precomputed product if its cardinality (max number of element of grade k) is lower than this threshold.
    // Else, the product is performed recursively (this is only for high dimensional Geometric Algebra).
    unsigned int maxDimPrecomputedProducts;

    // Each multivector component of grade k can have a dedicated basis accessor constant (E123 in "mv[E123]=42;") if its cardinality (max number of element of grade k) is lower than this threshold.
    // Else, this accessor is not created (this is only for high dimensional Geometric Algebra).
    unsigned int maxDimBasisAccessor;

    // numerical threshold for the numerical clean up
    double epsilon;

    // diagonal form of the metric after eigen vector/value decomposition
    Eigen::VectorXd diagonalMetric;

    // eigen vectors of the decomposition of the metric
    Eigen::MatrixXd transformationMatrix;

    // eigen vectors of the decomposition of the metric
    Eigen::MatrixXd inverseTransformationMatrix;



public:

    MetaData();

    MetaData(const std::string &filename);

    ~MetaData();

    void display() const;

    bool checkConsistency() const;

    bool metricDiagonalization();
};


#endif // GARAGEN_METADATA_HPP__
