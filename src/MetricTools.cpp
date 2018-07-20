// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// MetricTools.hpp
// This file is part of the Garamon Generator.
// Authors: Stephane Breuils and Vincent Nozick
// Conctact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program


#include "MetricTools.hpp"


#include <Eigen/Eigenvalues>
#include <cmath>
#include <iostream>


bool isMatrixDiagonal(const Eigen::MatrixXd &A, const double epsilon){

    for(unsigned int i=0; i<(unsigned int)A.rows(); ++i)
        for(unsigned int j=0; j<(unsigned int)A.cols(); ++j) {
            if(i==j) continue;
            if(fabs(A(i,j)) > epsilon)
                return false;
        }
    return true;
}


bool isMatrixIdentity(const Eigen::MatrixXd &A, const double epsilon){

    if(A.rows() != A.cols())
        return false;

    for(unsigned int i=0; i<(unsigned int)A.rows(); ++i)
        for(unsigned int j=0; j<(unsigned int)A.cols(); ++j) {
            if( (i==j) && (fabs(A(i,j)-1) > epsilon) )
                return false;

            if( (i!=j) && (fabs(A(i,j)) > epsilon) )
                return false;
        }

    return true;
}


bool isMatrixPermutationOfDiagonal(const Eigen::MatrixXd &metric, const double epsilon){

 return isMatrixDiagonal(metric.transpose()*metric,epsilon);
}


unsigned int getRank(const Eigen::MatrixXd &metric){

    Eigen::FullPivLU<Eigen::MatrixXd> lu(metric);
    return lu.rank();
}


void eigenDecomposition(const Eigen::MatrixXd &M, Eigen::MatrixXd &P, Eigen::MatrixXd &A){

    Eigen::EigenSolver<Eigen::MatrixXd> eigensolv;
    eigensolv.compute(M,true);
    P = eigensolv.eigenvectors().real();
    A = Eigen::MatrixXd::Zero(M.rows(), M.cols());
    A.diagonal() = eigensolv.eigenvalues().real();
}


double minAbsNonZeroValue(Eigen::VectorXd x){

    // deal only with positive values
    x = x.cwiseAbs();

    // for first min value
    bool firstValue = true;
    double minVal = 0.0; //std::numeric_limits<double>::infinity();
    for(unsigned int i=0; i<(unsigned int)x.size(); ++i){

        // zero values should be ignored
        if(x(i) < std::numeric_limits<double>::epsilon())
            continue;

        // if it is the first non zero value
        if(firstValue){
            minVal = x(i);
            firstValue = false;
        }else{
            // looking for the lowest value
            if(x(i) < minVal)
                minVal = x(i);      
        }   
    }

    return minVal; // negative if only zero values in the vector
}


// put the per column smallest element of P equal to 1 and update A consequently
// if P initially contains many sqrt(0.5), they will (probably) be transform to 1
// at the end, P is not an othrogonal matrix anymore.
Eigen::MatrixXd eigenRefinement(Eigen::MatrixXd &P, Eigen::MatrixXd &D, Eigen::MatrixXd &Pinv){

    // initialize the scale matrix to Identity
    Eigen::MatrixXd scaleMatrix = Eigen::MatrixXd::Identity(D.rows(),D.cols());

    // when possible, put integers in P, update Pinv consequently
    for(int i=0; i<P.rows(); ++i){

        // convert the per line smallest non-zero value to 1
        double minVal = minAbsNonZeroValue(P.col(i));
        P.col(i)    /= minVal;
        Pinv.row(i) *= minVal;
        //D(i,i) = 1.0 / (minVal*minVal);
        scaleMatrix(i,i) = 1.0/minVal;
    }
    return scaleMatrix;
}


// check if the clean up of the matrices still leads to a set of matrices whose product is equal to the initial metric.
bool checkNumericalCleanUp(const Eigen::MatrixXd &M, const Eigen::MatrixXd &P, const Eigen::MatrixXd & A, const double epsilon){

    Eigen::MatrixXd M2 = P * A * P.transpose();
    Eigen::MatrixXd nullMatrix = (M - M2).cwiseAbs();

    for(unsigned int i=0; i<(unsigned int)M.rows(); ++i)
        for(unsigned int j=0; j<(unsigned int)M.cols(); ++j)
            if(nullMatrix(i,j) > epsilon)
                return false;

    return true;
}


// check if the clean up of the matrices still leads to a set of matrices whose product is equal to the initial metric.
bool checkNumericalCleanUp(const Eigen::MatrixXd &M,
                           const Eigen::MatrixXd &P,
                           const Eigen::MatrixXd &A,
                           const Eigen::MatrixXd &Pinv,
                           const double epsilon){

    Eigen::MatrixXd M2 = P * A * Pinv;
    Eigen::MatrixXd identity = (M - M2).cwiseAbs();

    for(unsigned int i=0; i<(unsigned int)M.rows(); ++i)
        for(unsigned int j=0; j<(unsigned int)M.cols(); ++j)
            if(identity(i,j) > epsilon)
                return false;

    return true;
}


// when pertinent, replace a supposed integer value by the nearby int, etc.
Eigen::SparseMatrix<double> numericalCleanUpSparse(const Eigen::MatrixXd &M, const double epsilon){

    Eigen::SparseMatrix<double> N(M.rows(),M.cols());

    for(unsigned int i=0; i<(unsigned int)M.rows(); ++i)
        for(unsigned int j=0; j<(unsigned int)M.cols(); ++j){

            // ignore the zero
            if(fabs(M(i,j)) < epsilon)
                continue;

            // round near integers to integers
            int val = std::lround(M(i,j));
            if( fabs(val - M(i,j)) < epsilon ){
                N.insert(i,j) = val;
                continue;
            }

            // if failed to round with integer
            // round with decimal
            for(double d=0; d<=10; ++d){
                double decimal = -0.5 + d/10.0;
                if( fabs(val + decimal - M(i,j)) < epsilon ){
                    N.insert(i,j) = val + decimal;
                    continue;
                }
            }
        }

    return N;
}


// for vectors: when pertinent, replace a supposed integer value by the nearby int, etc.
Eigen::VectorXd vectorNumericalCleanUp(const Eigen::VectorXd& original, const double epsilon){
    Eigen::VectorXd outputVector(original);
    for(unsigned int j=0; j<(unsigned int)outputVector.size(); ++j) {

        // ignore the zero
        if (fabs(original(j)) < epsilon) {
            outputVector(j) = 0.0;
            continue;
        }

        // round near integers to integers
        int val = std::lround(original(j));
        if (fabs(val - original(j)) < epsilon) {
            outputVector(j) = val;
            continue;
        }

        // some negative power of 2
        const int maxNegPower = 7; // 2^{-6} = 0.015625 or the algebra dimension
        const double step = pow(2,-maxNegPower);
        for(double x=val-0.5; x<=val+0.5; x+=step)
            if(fabs(x - original(j)) < epsilon) {
                outputVector(j) = x;
                continue;
            }

        // if failed to round with integer
        // round with decimal
        for(int d=0; d<=10; ++d) {
            double decimal = -0.5 + d / 10.0;
            if (fabs(val + decimal - original(j)) < epsilon) {
                outputVector(j) = val + decimal;
                continue;
            }
        }
    }
    return outputVector;
}


// for matrices: when pertinent, replace a supposed integer value by the nearby int, etc.
Eigen::MatrixXd numericalCleanUp(const Eigen::MatrixXd &M, const double epsilon){

    Eigen::MatrixXd N = M; //Eigen::MatrixXd::Zero(M.rows(),M.cols());

    for(unsigned int i=0; i<(unsigned int)M.rows(); ++i)
        N.row(i) = vectorNumericalCleanUp(N.row(i), epsilon);

    return N;
}



Eigen::SparseMatrix<double, Eigen::ColMajor> computePerGradeTransformationMatrix(const Eigen::MatrixXd &vectorTransformationMatrix,
                                                                                 const unsigned int dimension,
                                                                                 const unsigned int grade,
                                                                                 const double epsilon){
    // Generate grade-vectors
    // i.e. dimension=3
    // generate: (e1+e2+e3) (e1+e2+e3) (e1+e2+e3)
    std::vector<std::vector<unsigned int> > listOfVectors;
    for(unsigned int j=0;j<dimension;++j){
        std::vector<unsigned int> oneVector;
        for(unsigned int i=0; i<dimension; ++i)
            oneVector.push_back((unsigned int) 1 << i);
        listOfVectors.push_back(oneVector);
    }


    // Compute the set of combination whose length is grade, in dimension space (i.e. for grade 2, dim 3: (1,2) (1,3) (2,3) )
    std::vector<std::vector<unsigned int> > sequence = generateCombinations(dimension, grade);

    // init the resulting transformation matrix
    Eigen::SparseMatrix<double, Eigen::ColMajor> resultSparse(sequence.size(), sequence.size());
//    Eigen::MatrixXd kVectorMetric = Eigen::MatrixXd::Zero(sequence.size(), sequence.size()); // to remove


    // link the indices used for the XOR wedge with the index in the sequence of combinations using an array
    // i.e. for grade 2, dim 3 :                (1,2) (1,3) (2,3)
    // order induced by the xor:                scal, 1, 2, 12, 3, 13, 23, 123
    // return their position in the sequence:      -, -, -,  0, -,  1,  2, -
    std::vector<unsigned int> xorIndex2Combination = getSetOfCombinationsFromXorIndexation(dimension, sequence);

    // generate the transformation matrix
    for(unsigned int l=0; l<sequence.size(); ++l) { // for each element of grade "grade'

        std::vector<unsigned int> currentVector = listOfVectors[0];
        std::vector<double>       currentCoeffs; // first line of the transformation matrix

        // fill the currentCoeffs using the transformation
        for(int idxTMat=0; idxTMat<vectorTransformationMatrix.cols(); ++idxTMat)
            currentCoeffs.push_back(vectorTransformationMatrix(sequence[l][0],idxTMat));


        std::vector<unsigned int> nextVector;
        std::vector<unsigned int> result;

        for(unsigned int j=0;j<sequence[l].size()-1;++j) {
            // Contains the number of wedge to be computed,
            // (v1^v2)^v3 for example

            nextVector = listOfVectors[j+1];
            std::vector<unsigned int> resultIdxTmp; // result = mv1 ^ mvss2
            std::vector<double> resultCoeffTmp;     // coefficient of the result

            // wedge between the last k-vector and a vector to make a (k+1)-vector up to 'grade'
            for(unsigned int i=0; i<currentVector.size(); ++i) {
                for(unsigned int k=0; k<nextVector.size(); ++k) {

                    // index initialisation
                    unsigned int mvC = 0;
                    unsigned int mvA = currentVector[i];
                    unsigned int mvB = nextVector[k];

                    // coefficient computation
                    double coeffA = currentCoeffs[i];
                    double coeffB = vectorTransformationMatrix(sequence[l][j+1],k);
                    double coeffC = 0.0; // To be changed with the transformation matrix
                    if((coeffA != 0.0) && (coeffB != 0.0)){
                        getSignAndBladeOuterProduct(mvC, mvA, mvB, coeffC, coeffA, coeffB);
                        // update the results
                        if(coeffC !=0.0) {
                            resultIdxTmp.push_back(mvC);
                            resultCoeffTmp.push_back(coeffC);
                        }
                    }
                }
            }

            currentVector = resultIdxTmp;
            currentCoeffs = resultCoeffTmp;
            resultIdxTmp.clear();
            resultCoeffTmp.clear();
        }

        for(unsigned int idxRes=0; idxRes<currentVector.size(); ++idxRes){
            if(fabs(currentCoeffs[idxRes])>epsilon)
                resultSparse.coeffRef(l,xorIndex2Combination[currentVector[idxRes]]) += currentCoeffs[idxRes];
        }

    }
    return resultSparse;
}


// convert the grade transformation matrix (inverse or not) to a vector of values (row, colums,value). no interpretation are done
std::vector<double> transformationMatricesToVectorOfComponents(
        const Eigen::SparseMatrix<double, Eigen::ColMajor> &transformationMatrix, const int grade, const bool isInverse) {
    std::vector<double> outputVector;
    for (int k = 0; k < transformationMatrix.outerSize(); ++k)
        for (Eigen::SparseMatrix<double, Eigen::ColMajor>::InnerIterator it(transformationMatrix,k); it; ++it) {
            outputVector.push_back((double)it.row());
            outputVector.push_back((double)it.col());
            outputVector.push_back((double)it.value());
        }
    return outputVector;
}


// Compute the inverse of a transformation matrix given by transformationMatrix
Eigen::SparseMatrix<double, Eigen::ColMajor>  computeInverseTransformationMatrix(const Eigen::SparseMatrix<double, Eigen::ColMajor>& transformationMatrix, const double epsilon){
    Eigen::SparseMatrix<double, Eigen::ColMajor> inverseTransformationMatrix(transformationMatrix.rows(),transformationMatrix.rows());
    Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>> solver;
    solver.analyzePattern(transformationMatrix);
    solver.factorize(transformationMatrix);
    // compute the inverse M.M-1 = Id
    // i.e. solve all system where the unknown vector is a column of M-1
    for(unsigned int j=0; j<(unsigned int)inverseTransformationMatrix.cols(); ++j){
        Eigen::SparseVector<double> vecId(inverseTransformationMatrix.rows());
        vecId.insert(j) = 1.0;
        Eigen::VectorXd x = solver.solve(vecId);
        // recopy into sparse matrix
        for(unsigned i=0; i<(unsigned int)x.size(); ++i)
            if(fabs(x(i)) > epsilon)
                inverseTransformationMatrix.insert(i,j) = x(i);
    }

    return inverseTransformationMatrix;
}


// Compute the transformation matrices for grades ranging from 1 to d
// For each per-grade transformation matrix, we construct a string containing the list of its non-zero elements.
// We also pick up the number of non zero elements of each sparse matrices, the result is put into transformationMatricesSize
std::pair<std::vector<double>,std::vector<double>> computeTransformationMatricesToVector(const Eigen::MatrixXd &P, const double epsilon, std::vector<unsigned int>& transformationMatricesSizes,
                                                                                       std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor> >& allTransformationMatrices,
                                                                                       std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor> >& allInverseTransformationMatrices){
    std::pair<std::vector<double>,std::vector<double>> transformationMatrices; // contains non-inverse and inverse transformation matrices

    // push the scalar transformation matrix
    Eigen::SparseMatrix<double, Eigen::ColMajor> spscalarTransformationMatrix(1, 1);
    spscalarTransformationMatrix.insert(0,0)=1.0;
    allTransformationMatrices.push_back(spscalarTransformationMatrix);

    // The number of non-zeros elements of the grade 0 transformation matrix is simply 1
    transformationMatricesSizes.push_back(1);

    // SCALAR transformation matrix to a list of triplets (here only 1) which will be generated afterwards
    std::vector<double> currentBasisTransformComponents = transformationMatricesToVectorOfComponents(spscalarTransformationMatrix, 0, false);
    transformationMatrices.first.insert(std::end(transformationMatrices.first), std::begin(currentBasisTransformComponents), std::end(currentBasisTransformComponents));
    currentBasisTransformComponents = transformationMatricesToVectorOfComponents(spscalarTransformationMatrix, 0, true); // inverse transformation matrices
    transformationMatrices.second.insert(std::end(transformationMatrices.second), std::begin(currentBasisTransformComponents), std::end(currentBasisTransformComponents));
    allInverseTransformationMatrices.push_back(spscalarTransformationMatrix);

    // push the transformation matrix to the vector of transformation matrices
    Eigen::SparseMatrix<double, Eigen::ColMajor> spVectorTransformationMatrix(P.rows(), P.cols());
    for(unsigned int j=0; j<(unsigned int)P.cols(); ++j)
        for(unsigned int i=0; i<(unsigned int)P.rows(); ++i)
            if(fabs(P(i,j)) >epsilon) // if non-zero value, insert it
                spVectorTransformationMatrix.insert(i,j)=P(i,j);
    allTransformationMatrices.push_back(spVectorTransformationMatrix);

    // extract the number of non-zeros elements of the grade 1transformation matrix
    transformationMatricesSizes.push_back((unsigned int)spVectorTransformationMatrix.nonZeros());

    // convert the grade 1 transformation matrix to a list of triplets which will be generated afterwards
    currentBasisTransformComponents = transformationMatricesToVectorOfComponents(spVectorTransformationMatrix, 0, false);
    transformationMatrices.first.insert(std::end(transformationMatrices.first), std::begin(currentBasisTransformComponents), std::end(currentBasisTransformComponents));

    // compute also the inverse transformation matrix
    Eigen::SparseMatrix<double, Eigen::ColMajor> spVectorInverseTransformation = computeInverseTransformationMatrix(spVectorTransformationMatrix, epsilon);

    // and convert it to a list of triplets
    currentBasisTransformComponents = transformationMatricesToVectorOfComponents(spVectorInverseTransformation, 0, true);
    transformationMatrices.second.insert(std::end(transformationMatrices.second), std::begin(currentBasisTransformComponents), std::end(currentBasisTransformComponents));
    allInverseTransformationMatrices.push_back(spVectorInverseTransformation);

    // remaining transformation matrices
    for(unsigned int i=2;i<=(unsigned int)P.cols();++i){

        Eigen::SparseMatrix<double, Eigen::ColMajor> spPerGradeTransformation = computePerGradeTransformationMatrix(P, (unsigned int) P.rows(), i, epsilon); // resultMatrix, inputMatrix, dimension, nb vectors to wedge

        // extract the number of non-zeros elements of the current transformation matrix
        transformationMatricesSizes.push_back((unsigned int)spPerGradeTransformation.nonZeros());

        // convert the current transformation matrix to a list of triplets which will be generated afterwards
        currentBasisTransformComponents = transformationMatricesToVectorOfComponents(spPerGradeTransformation, 0, false);
        transformationMatrices.first.insert(std::end(transformationMatrices.first), std::begin(currentBasisTransformComponents), std::end(currentBasisTransformComponents));
        allTransformationMatrices.push_back(spPerGradeTransformation);

        // compute also the inverse of the current transformation matrix
        Eigen::SparseMatrix<double, Eigen::ColMajor> spPerGradeInverseTransformation = computeInverseTransformationMatrix(spPerGradeTransformation, epsilon);

        // convert the current transformation matrix to a list of triplets which will be generated afterwards
        currentBasisTransformComponents = transformationMatricesToVectorOfComponents(spPerGradeInverseTransformation, 0, true);
        transformationMatrices.second.insert(std::end(transformationMatrices.second), std::begin(currentBasisTransformComponents), std::end(currentBasisTransformComponents));
        allInverseTransformationMatrices.push_back(spPerGradeInverseTransformation);
    }
    return transformationMatrices;
}
