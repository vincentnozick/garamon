// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// MetaData.cpp
// This file is part of the Garamon Generator.
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program


#include <algorithm>
#include <cctype>
#include <cmath>
#include <iostream>
#include <ostream>
#include <set>
#include <stdexcept>
#include <string>

#include <Eigen/Core>

#include "MetaData.hpp"
#include "ConfigParser.hpp"
#include "MetricTools.hpp"


MetaData::MetaData() {}

MetaData::~MetaData() {}

void MetaData::display() const {
    std::cout << "MetaData" << std::endl;
    std::cout << "dimension         : " << dimension << std::endl;
    std::cout << "namespace         : " << namespaceName << std::endl;

    std::cout << "refinement        : " << (useEigenRefinement==true?"true":"false") << std::endl;
    std::cout << "cleanup           : " << (useNumericalCleanUp==true?"true":"false") << std::endl;
    std::cout << "maxDim prec func  : " << maxDimPrecomputedProducts << std::endl;
    std::cout << "maxDim accessors  : " << maxDimBasisAccessor << std::endl;
    std::cout << "epsilon           : " << epsilon << std::endl;
    std::cout << "basis vector name : ";
    for(unsigned int i=0; i<basisVectorName.size(); ++i)
        std::cout << basisVectorName[i] << " ";
    std::cout << std::endl;
    std::cout << "metric            : \n" << metric << std::endl;
    std::cout << "is initialy diag  : " << (inputMetricDiagonal==true?"true":"false") << std::endl;
    std::cout << "is identity       : " << (identityMetric==true?"true":"false") << std::endl;
    std::cout << "is full rank      : " << (fullRankMetric==true?"true":"false") << std::endl;
    std::cout << "is permutation of diagonal matrix : " << (inputMetricPermutationOfDiagonal==true?"true":"false") << std::endl;
    std::cout << "diag metric       : " << diagonalMetric.transpose() << std::endl;
    if(!inputMetricDiagonal){
        std::cout << "Vector transformation matrix : \n" << transformationMatrix << std::endl;
        std::cout << "Vector inverse Transformation  : \n" << inverseTransformationMatrix<< std::endl;
    }
}

bool MetaData::checkConsistency() const {

    bool consistencyCheck = true;

    // maxDimBasisAccessor consistency: at least for vectors
    if(maxDimBasisAccessor == 0){
        std::cout << "error: at least vector accessors are required, see 'max dimension basis accessor' in the conf file." << std::endl;
        consistencyCheck = false;
    }

    // dimension consistency
    if(dimension == 0){
        std::cout << "error: dimension should not be 0." << std::endl;
        consistencyCheck = false;
    }

    // the multivectors are indexed with 32 bits xor indices
    if(dimension > 31){
        std::cout << "error: dimension should not be higher than 31." << std::endl;
        consistencyCheck = false;
    }

    // basis vector name dimension consistency
    bool basisNamesConsistent = true;
    if(basisVectorName.size() != dimension){
        std::cout << "error: 'basis vector name' size is not consistent with 'dimension'." << std::endl;
        consistencyCheck = false;
        basisNamesConsistent = false;
    }

    // metric empty
    bool metricDefined = true;
    if( (metric.rows()==0) || (metric.cols()==0) ) {
        metricDefined = false;
        std::cout << "error: the metric matrix is not defined." << std::endl;
        consistencyCheck = false;
    }

    // metric matrix ?
    if(metricDefined) {
        // metric : square matrix ?
        bool squareMatrix = true;
        if( metric.rows() != metric.cols() ){
            std::cout << "error: the metric is not a square matrix." << std::endl;
            consistencyCheck = false;
            squareMatrix = false;
        }

        // metric : dimension consistency
        if( ((unsigned int)metric.cols() != dimension) || ((unsigned int)metric.rows() != dimension) ) {
            std::cout << "error: the metric dimension is not consistent with 'dimension'." << std::endl;
            consistencyCheck = false;
        }

        if(squareMatrix) {
            // check if metric is symetric
            bool symetric = true;
            for(unsigned int i=0; i<(unsigned int)metric.rows(); ++i)
                for(unsigned int j=i; j<(unsigned int)metric.cols(); ++j)
                    if(std::fabs(metric(i,j) - metric(j,i)) > epsilon)
                        symetric = false;
            if (!symetric) {
                std::cout << "error: the metric is not a symmetric matrix." << std::endl;
                consistencyCheck = false;
            }
        }
    }

    // check if the name of the vector basis are not ambiguous for high dimensions
    // i.e. in dimension 15: e12 is for twelve or one-two ?
    // (only possible if there is one name per basis vector)
    if(basisNamesConsistent && dimension > 0 && dimension <= 31) {
        std::set<std::string> basis;

        // represent a k-vector with a binary number (k-st bit to 1 means the k-st basis is used)
        for(unsigned int i=1; i < (1u << dimension); ++i){

            std::string kvector;
            for(unsigned int k=0; k<dimension; ++k)
                if(i & (1u << k))
                    kvector = kvector + basisVectorName[k];

            if(basis.count(kvector) != 0){
                std::cout << "error in the basis vector name: " << kvector << " is ambiguous." << std::endl;
                consistencyCheck = false;
            }else{
                basis.insert(kvector);
            }
        }

        // the basis vector names are part of C++ identifiers (e.g. E12, e12())
        for(const auto &name : basisVectorName)
            if(!std::all_of(name.begin(), name.end(), [](unsigned char c){return std::isalnum(c) || c == '_';})){
                std::cout << "error in the basis vector name: '" << name << "' should contain only letters, digits or '_'." << std::endl;
                consistencyCheck = false;
            }
    }

    // namespace name compatible with C++
    if(namespaceName.empty() || std::isalpha((unsigned char)namespaceName[0]) == 0){
        std::cout << "error in the namespace name: the first character of '" << namespaceName << "' should be an alphabetic letter for C++ compliance." << std::endl;
        consistencyCheck = false;
    }else if(!std::all_of(namespaceName.begin(), namespaceName.end(), [](unsigned char c){return std::isalnum(c) || c == '_';})){
        std::cout << "error in the namespace name: '" << namespaceName << "' should contain only letters, digits or '_' for C++ compliance." << std::endl;
        consistencyCheck = false;
    }

    return consistencyCheck;
}

MetaData::MetaData(const std::string &filename) : MetaData(ConfigParser(filename)) {}

MetaData::MetaData(const ConfigParser &parser) {

    // load all components of the meta data
    if(!parser.readString("namespace", namespaceName))
        throw std::runtime_error("failed to find namespace");

    if(!parser.readUInt("dimension", dimension))
        throw std::runtime_error("failed to find dimension");

    if(!parser.readUInt("max dimension precomputed products", maxDimPrecomputedProducts))
        throw std::runtime_error("failed to find max dimension precomputed products");

    if(!parser.readUInt("max dimension basis accessor", maxDimBasisAccessor))
        throw std::runtime_error("failed to find max dimension basis accessor");

    if(!parser.readStringList("basis vector name", basisVectorName))
        throw std::runtime_error("failed to find basis vector name");

    if(!parser.readMatrix("metric", metric))
        throw std::runtime_error("failed to find metric");

    if(!parser.readBool("metric decomposition refinement", useEigenRefinement))
        throw std::runtime_error("failed to find metric decomposition refinement");

    if(!parser.readBool("metric decomposition numerical cleanup", useNumericalCleanUp))
        throw std::runtime_error("failed to find metric decomposition numerical cleanup");

    if(!parser.readDouble("metric decomposition numerical cleanup espilon", epsilon))
        throw std::runtime_error("failed to find metric decomposition numerical cleanup espilon");

    // check the data consistency
    if(!checkConsistency())
        throw std::runtime_error("configuration inconsistent ... abort");

    // metric diagonalization
    if(!metricDiagonalization())
        throw std::runtime_error("metric diagonalization ... failed\n"
                                 "Please try again without the metric decomposition refinement or without the numerical cleanup.");
}


bool MetaData::metricDiagonalization() {

	// compute the metric rank
	if (getRank(metric) == dimension)
		fullRankMetric = true;
	else fullRankMetric = false;

    // check if the metric is already diagonal
    if(isMatrixDiagonal(metric, epsilon)){
        diagonalMetric = metric.diagonal();
        inputMetricDiagonal = true;
        if(isMatrixIdentity(metric,epsilon))
            identityMetric = true;
        return true;
    }


    // ckeck if the metric is a permutation of a diagonal matrix (for fast dual)
    inputMetricPermutationOfDiagonal = isMatrixPermutationOfDiagonal(metric, epsilon);

    // compute the diagonalization
    Eigen::MatrixXd diagonalMatrix;
    eigenDecomposition(metric, transformationMatrix, diagonalMatrix);
//    std::cout << "metric\n" << metric << std::endl;
//    std::cout << "P\n" << transformationMatrix << std::endl;
//    std::cout << "D\n" << diagonalMatrix << std::endl;
//    std::cout << "metric ??\n" << transformationMatrix * diagonalMatrix * transformationMatrix.transpose() << std::endl;

    // invert the transformation matrix
    inverseTransformationMatrix = Eigen::MatrixXd(transformationMatrix.transpose());

    // convert floating points to nearest integers (when possible)
    Eigen::MatrixXd scaleMatrix = Eigen::MatrixXd::Identity(metric.rows(),metric.cols()); // refer to issue #9
    if(useEigenRefinement)
        scaleMatrix = eigenRefinement(transformationMatrix, diagonalMatrix, inverseTransformationMatrix);
//    std::cout << "metric\n" << metric << std::endl;
//    std::cout << "P\n" << transformationMatrix << std::endl;
//    std::cout << "D\n" << diagonalMatrix << std::endl;
//    std::cout << "Pinv\n" << inverseTransformationMatrix << std::endl;
//    std::cout << "metric ??\n" << transformationMatrix * diagonalMatrix * inverseTransformationMatrix << std::endl;


    //numerical clean up
    if(useNumericalCleanUp) {
        transformationMatrix = numericalCleanUp(transformationMatrix,epsilon);
        diagonalMatrix = numericalCleanUp(diagonalMatrix, epsilon);
        inverseTransformationMatrix = numericalCleanUp(inverseTransformationMatrix,epsilon);
    }

    // check if the new metric is identity
    if(isMatrixIdentity(diagonalMatrix,epsilon))
        identityMetric = true;

    // check decomposition
    if(!checkNumericalCleanUp(metric,transformationMatrix,diagonalMatrix,inverseTransformationMatrix, epsilon))
        return false;

    // diagonal metric must be changed if we consider only inverse
    diagonalMatrix = (scaleMatrix*scaleMatrix)*diagonalMatrix;
    if(useNumericalCleanUp)
        diagonalMatrix = numericalCleanUp(diagonalMatrix, epsilon);

    // compact form of the diagonal matrix
    diagonalMetric = diagonalMatrix.diagonal();

    return true;
}

