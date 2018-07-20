// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// MetaData.cpp
// This file is part of the Garamon Generator.
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program


#include <math.h>
#include <set>

#include <Eigen/Eigenvalues>

#include "MetaData.hpp"
#include "ConfigParser.hpp"


MetaData::MetaData() : dimension(0), inputMetricDiagonal(false), identityMetric(false), inputMetricPermutationOfDiagonal(false), maxDimPrecomputedProducts(256) {}

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

    // basis vector name dimension consistency
    if(basisVectorName.size() != dimension){
        std::cout << "error: 'basis vector name' size is not consistent with 'dimension'." << std::endl;
        consistencyCheck = false;
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
                    if(fabs(metric(i,j) - metric(j,i)) > epsilon)
                        symetric = false;
            if (!symetric) {
                std::cout << "error: the metric is not a symmetric matrix." << std::endl;
                consistencyCheck = false;
            }
        }
    }

    // check if the name of the vector basis are not ambiguous for high dimensions
    // i.e. in dimension 15: e12 is for twelve or one-two ?
    std::set<std::string> basis;

    // represent a k-vector with a binary number (k-st bit to 1 means the k-st basis is used)
    for(unsigned int i=1; i<=pow(2,dimension); ++i){

        std::string kvector;
        for(unsigned int k=0; k<dimension; ++k)
            if(i & (1 << k))
                kvector = kvector + basisVectorName[k];

        if(basis.count(kvector) != 0){
            std::cout << "error in the basis vector name: " << kvector << " is ambiguous." << std::endl;
            consistencyCheck = false;
        }else{
            basis.insert(kvector);
        }
    }

    // namespace name compatible with C++
    if(isalpha(namespaceName[0]) == 0){
        std::cout << "error in the namespace name: the first character of '" << namespaceName << "' should be an alphabetic letter for C++ compliance." << std::endl;
        consistencyCheck = false;
    }

    return consistencyCheck;
}

MetaData::MetaData(const std::string &filename):inputMetricDiagonal(false), identityMetric(false) {

    // open the parser
    ConfigParser parser(filename);

    // load all components of the meta data
    if(!parser.readString("namespace", namespaceName)) {
        std::cerr << "error: failed to find " << "namespace" << std::endl;
        exit(EXIT_FAILURE);
    }

    if(!parser.readUInt("dimension", dimension)){
        std::cerr << "error: failed to find " << "dimension" << std::endl;
        exit(EXIT_FAILURE);
    }

    if(!parser.readUInt("max dimension precomputed products", maxDimPrecomputedProducts)){
        std::cerr << "error: failed to find " << "max dimension precomputed products" << std::endl;
        exit(EXIT_FAILURE);
    }

    if(!parser.readUInt("max dimension basis accessor", maxDimBasisAccessor)){
        std::cerr << "error: failed to find " << "max dimension basis accessor" << std::endl;
        exit(EXIT_FAILURE);
    }

    if(!parser.readStringList("basis vector name", basisVectorName)){
        std::cerr << "error: failed to find " << "basis vector name" << std::endl;
        exit(EXIT_FAILURE);
    }

    if(!parser.readMatrix("metric", metric)){
        std::cerr << "error: failed to find " << "metric" << std::endl;
        exit(EXIT_FAILURE);
    }

    if(!parser.readBool("metric decomposition refinement", useEigenRefinement)){
        std::cerr << "error: failed to find " << "metric decomposition refinement" << std::endl;
        exit(EXIT_FAILURE);
    }

    if(!parser.readBool("metric decomposition numerical cleanup", useNumericalCleanUp)){
        std::cerr << "error: failed to find " << "metric decomposition numerical cleanup" << std::endl;
        exit(EXIT_FAILURE);
    }

    if(!parser.readDouble("metric decomposition numerical cleanup espilon", epsilon)){
        std::cerr << "error: failed to find " << "metric decomposition numerical cleanup espilon" << std::endl;
        exit(EXIT_FAILURE);
    }

    // check the data consistency
    if(!checkConsistency()){
        std::cerr << "configuration inconsistent ... abord" << std::endl;
        exit(EXIT_FAILURE);
    }

    // metric diagonalization
    if(!metricDiagonalization()){
        std::cerr << "metric diagonalization ... failed" << std::endl;
        std::cerr << "Please try again without the metric decomposition refinement or without the numerical cleanup." << std::endl;
        exit(EXIT_FAILURE);
    }
}


bool MetaData::metricDiagonalization() {

    // check if the metric is already diagonal
    if(isMatrixDiagonal(metric, epsilon)){
        diagonalMetric = metric.diagonal();
        inputMetricDiagonal = true;
        if(isMatrixIdentity(metric,epsilon))
            identityMetric = true;
        return true;
    }

    // compute the metric rank
    if(getRank(metric) == dimension)
        fullRankMetric = true;
    else fullRankMetric = false;

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
    Eigen::MatrixXd scaleMatrix(metric.rows(),metric.cols());
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

