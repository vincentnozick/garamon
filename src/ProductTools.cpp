// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// ProductTools.cpp
// This file is part of the Garamon Generator.
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

#include "ProductTools.hpp"

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdlib>
#include <list>
#include <map>
#include <vector>

#include <Eigen/Core>
#include <Eigen/LU> // IWYU pragma: keep (MatrixBase::determinant() is declared in Core but defined in LU)
#include <Eigen/SparseCore>


// fills the [xor to pos/grade] and the [pos/grade to xor] arrays.
ProductTools::ProductTools(const unsigned int vectorSpaceDimension)
        : dimension(vectorSpaceDimension), xorIndexToGradeAndPosition(std::size_t(1) << vectorSpaceDimension)
{
    // fill simultaneously both 'gradePositionToXorIndex' and 'xorIndexToGradeAndPosition'
    for(unsigned int grade=0; grade<=dimension; ++grade) {

        // for the current grade, generate all the XOR indices
        std::vector<bool> booleanXorIndex(dimension);
        std::fill(booleanXorIndex.begin(), booleanXorIndex.end(), false);

        // build the first xorIndex of grade 'grade' (with 'grade' values to true).
        std::fill(booleanXorIndex.begin(), booleanXorIndex.begin() + grade, true);

        // initialize a new grade entry for 'gradePositionToXorIndex'
        gradePositionToXorIndex.push_back(std::vector<unsigned int>());

        // for each permutation of the true values of booleanXorIndex
        unsigned int positionInKVector = 0;
        do {

            // convert the vector of bool to a unsigned int (each bit get a boolean value)
            unsigned int xorIndex = 0;
            for(unsigned int i=0; i<dimension; ++i) {
                if(booleanXorIndex[i]) {
                    xorIndex += (1<<i); // Xor index
                }
            }

            // update 'gradePositionToXorIndex' array (for the current grade, for the last pos)
            gradePositionToXorIndex.back().push_back(xorIndex);

            // update 'xorIndexToGradeAndPosition' array
            xorIndexToGradeAndPosition.at(xorIndex) = {grade,positionInKVector};
            positionInKVector++;

        } while(std::prev_permutation(booleanXorIndex.begin(), booleanXorIndex.end()));
    }
}


ProductTools::~ProductTools() {}


unsigned int ProductTools::hammingWeight(const unsigned int xorIndexMv) {
    return (unsigned int)std::popcount(xorIndexMv);
}


bool ProductTools::outerProductExists(const unsigned int xorIndexMv1, const unsigned int xorIndexMv2) {
    return (xorIndexMv1 & xorIndexMv2) == 0;
}


int ProductTools::outerProductSign(const unsigned int xorIndexMv1, const unsigned int xorIndexMv2)
{
    unsigned int nbOfOnes = 0;
    unsigned int tmpMv = xorIndexMv1;
    tmpMv = tmpMv >> 1;
    while(tmpMv!=0){
        // Count number of ones in common
        nbOfOnes += hammingWeight(tmpMv & xorIndexMv2);
        tmpMv = tmpMv >> 1;
    }
    return (-2*(nbOfOnes%2))+1; // (-1)^{nbOfOnes}
}


bool ProductTools::innerProductExists(const unsigned int xorIndexMv1, const unsigned int xorIndexMv2)
{
    // compute the grade of the 2 multivectors
    unsigned int gradeMv1 = hammingWeight(xorIndexMv1);
    unsigned int gradeMv2 = hammingWeight(xorIndexMv2);

    // Left contraction contribution
    if((gradeMv1 < gradeMv2)){
        if((xorIndexMv1 & (~(xorIndexMv2))) == 0){
            // the left contraction part may exist
            return true;
        }
    // the right contraction contribution
    }else if (gradeMv1 > gradeMv2){
        if((xorIndexMv2 & (~(xorIndexMv1))) == 0){
            // the right contraction part may exist
            return true;
        }
    // scalar product contribution
    }else{
        // gradeMv1 == gradeMv2
        if((xorIndexMv1 ^ xorIndexMv2) == 0) {
            // the scalar product part may exist
            return true;
        }
    }

    // this inner product does not lead to any computation / results
    return false;
}



unsigned int ProductTools::productResultXorIndex(const unsigned int xorIndexMv1, const unsigned int xorIndexMv2){
    return xorIndexMv1 ^ xorIndexMv2;
}


// The inner and geometric product are a metric dependant so we have to take the metric into account.
double ProductTools::productCoefficientFromMetric(const unsigned int xorIndexMv1, const unsigned int xorIndexMv2,
                                                  const Eigen::VectorXd &diagonalMetric){

    // First determine the sign, ie. number of permutation needed to have the common representation of mv1 and mv2 ordered.
    double coefficient = outerProductSign(xorIndexMv1,xorIndexMv2);

    // Each common bit of xorIndexMv1 and xorIndexMv2 will imply the corresponding metric value.
    // Example: suppose we compute the inner product of e123 (xorIndex: 0111) with e23 (xorIndex 0110) in 4D.
    // The commons 1s is 0110 (corresponding to e23).
    // thus coefficient is multiplied by metric[0100]*metric[0010]
    unsigned int commonOnes = xorIndexMv1 & xorIndexMv2;

    // loop over all the common 1s bit of xorIndexMv1 and xorIndexMv2
    Eigen::Index dimension = diagonalMetric.size();
    for(Eigen::Index i=0;i<dimension; ++i){
        if((commonOnes & (1<<i)) != 0){
            // one at position i in the binary representation of commonOnes
            coefficient *= diagonalMetric[i];
        }
    }

    return coefficient;
}


std::list<productComponent<double>> ProductTools::generateExplicitOuterProductList(const unsigned int gradeMv1, const unsigned int gradeMv2) const{

    std::list<productComponent<double>> outputListOfProducts;

    // for each element of grade 'gradeMv1'
    for(auto indexMv1 : gradePositionToXorIndex[gradeMv1])
        // for each element of grade 'gradeMv2'
        for(auto indexMv2 : gradePositionToXorIndex[gradeMv2]){

            // Compute outer product of the blades represented by the index indexMv1 and indexMv2
            productComponent<double> currentProduct;
            currentProduct.indexOfMv1 = xorIndexToGradeAndPosition[indexMv1].second;
            currentProduct.indexOfMv2 = xorIndexToGradeAndPosition[indexMv2].second;

            // if indexMv1 and indexMv2 has no common vector basis, we can compute their outer product
            if(outerProductExists(indexMv1,indexMv2)){

                // compute the product data
                currentProduct.coefficient = outerProductSign(indexMv1,indexMv2);
                unsigned int indexOfResult = productResultXorIndex(indexMv1, indexMv2); // compute indexMv1^indexMv2
                currentProduct.indexOfMv3  = xorIndexToGradeAndPosition[indexOfResult].second;

                // add this product to the fine component of mv3
                outputListOfProducts.push_back(currentProduct);
            }
        }

    return outputListOfProducts;
}


std::list<productComponent<double>> ProductTools::generateExplicitInnerProductListEuclideanSpace(const unsigned int gradeMv1, const unsigned int gradeMv2,
                                                                                            const Eigen::VectorXd &diagonalMetric) const{

    std::list<productComponent<double>> outputListOfProducts;

    // for each element of grade 'gradeMv1'
    for(auto indexMv1 : gradePositionToXorIndex[gradeMv1])
        // for each element of grade 'gradeMv2'
        for(auto indexMv2 : gradePositionToXorIndex[gradeMv2]){

            // Compute outer product of the blades represented by the index indexMv1 and indexMv2
            productComponent<double> currentProduct;
            currentProduct.indexOfMv1 = xorIndexToGradeAndPosition[indexMv1].second;
            currentProduct.indexOfMv2 = xorIndexToGradeAndPosition[indexMv2].second;

            // compute the inner product
            if(innerProductExists(indexMv1, indexMv2)){
                // compute the coefficient
                currentProduct.coefficient = productCoefficientFromMetric(indexMv1, indexMv2,diagonalMetric);

                // degenerate metric: the product may be null
                if(currentProduct.coefficient == 0.0)
                    continue;
                // compute the resulting index
                unsigned int indexOfResult = xorIndexToGradeAndPosition[productResultXorIndex(indexMv1, indexMv2)].second;
                // fill the resulting (gradeMv2-gradeMv1)-vector which is expressed in the orthogonal basis
                currentProduct.indexOfMv3  = indexOfResult;



                // add this product to the fine component of mv3
                outputListOfProducts.push_back(currentProduct);
            }
        }

    return outputListOfProducts;
}


std::list<productComponent<double>> ProductTools::generateExplicitInnerProductList(const unsigned int gradeMv1,
                                                                                   const unsigned int gradeMv2,
                                                                                   const Eigen::SparseMatrix<double, 0> &transformationMatrixMv1,
                                                                                   const Eigen::SparseMatrix<double, 0> &transformationMatrixMv2,
                                                                                   const Eigen::SparseMatrix<double, 0> &transformationMatrixMv3,
                                                                                   const Eigen::VectorXd &diagonalMetric) const{
    std::list<productComponent<double>> outputListOfProducts;

    // for each element of grade 'gradeMv1'
    for(unsigned int indexMv1=0; indexMv1<(unsigned int)transformationMatrixMv1.rows(); ++indexMv1){

        // for each element of grade 'gradeMv2'
        for(unsigned int indexMv2=0; indexMv2<(unsigned int)transformationMatrixMv2.rows(); ++indexMv2){

            // create the vector results according to its number of elements
            Eigen::SparseVector<double,Eigen::ColMajor> mv3Ortho(transformationMatrixMv3.cols());

            // compute the coefficient of the current product, in the orthogonal basis
            // -> for each basis vector in the orthogonal space of grade 'gradeMv1'
            for(Eigen::SparseMatrix<double>::InnerIterator itMv1(transformationMatrixMv1,indexMv1); itMv1; ++itMv1){

                // -> for each basis vector in the orthogonal space of grade 'gradeMv2'
                for(Eigen::SparseMatrix<double>::InnerIterator itMv2(transformationMatrixMv2,indexMv2); itMv2; ++itMv2){

                    // get the xor indices of Mv1 and Mv2 (independent of the basis)
                    unsigned int xorIndexMv1 = gradePositionToXorIndex[gradeMv1][itMv1.index()];
                    unsigned int xorIndexMv2 = gradePositionToXorIndex[gradeMv2][itMv2.index()];

                    // compute the inner product
                    if(innerProductExists(xorIndexMv1, xorIndexMv2)){
                        // compute the coefficient
                        double coefficient = productCoefficientFromMetric(xorIndexMv1, xorIndexMv2,diagonalMetric);

                        // compute the resulting index
                        unsigned int indexOfResult = xorIndexToGradeAndPosition[productResultXorIndex(xorIndexMv1, xorIndexMv2)].second;

                        // fill the resulting (gradeMv2-gradeMv1)-vector which is expressed in the orthogonal basis
                        double resultingCoefficient = coefficient * itMv1.value() * itMv2.value();

                        // in a sparse vector, add only non-zero values
                        if(resultingCoefficient !=0)
                            mv3Ortho.coeffRef(indexOfResult) += resultingCoefficient;
                    }
                }
            } // end of: compute the coefficient of the current product, in the orthogonal basis

            // convert back the result in the original basis
            Eigen::SparseVector<double,Eigen::ColMajor> mv3 = transformationMatrixMv3 * mv3Ortho;

            // add all non-zero component to the list of products
            for(Eigen::SparseVector<double>::InnerIterator itMv3(mv3); itMv3; ++itMv3) {
                productComponent<double> currentProduct;
                if(itMv3.value() != 0.0){      // maybe useless test, but we don't care ...
                    currentProduct.indexOfMv1  = indexMv1;
                    currentProduct.indexOfMv2  = indexMv2;
                    currentProduct.indexOfMv3  = itMv3.index();
                    currentProduct.coefficient = itMv3.value();
                    outputListOfProducts.push_back(currentProduct);
                }
            }
        }
    }
    return outputListOfProducts;
}


std::vector<std::list<productComponent<double>>> ProductTools::generateExplicitGeometricProductListEuclideanSpace(const unsigned int gradeMv1,const unsigned int gradeMv2, const Eigen::VectorXd& diagonalMetric) const{

    std::vector<std::list<productComponent<double>> > outputListOfProducts(dimension+1);
    // As we ignore the inner and outer products contributions. Thus, all products whose resulting grade is gradeOuter=gradeMv1+gradeMv2 or gradeInner=abs(gradeMv1-gradeMv2) will be ignored
    const unsigned int gradeOuter = gradeMv1+gradeMv2;
    const auto gradeInner = (unsigned int)std::abs((float)gradeMv1-gradeMv2);

    // for each element of grade 'gradeMv1'
    for(auto indexMv1 : gradePositionToXorIndex[gradeMv1])
        // for each element of grade 'gradeMv2'
        for(auto indexMv2 : gradePositionToXorIndex[gradeMv2]){

            // Compute outer product of the blades represented by the index indexMv1 and indexMv2
            productComponent<double> currentProduct;

                // compute the coefficient
                double coefficient = productCoefficientFromMetric(indexMv1, indexMv2,diagonalMetric);

                // compute the resulting index
                unsigned int indexOfResult = xorIndexToGradeAndPosition[productResultXorIndex(indexMv1, indexMv2)].second;
                // compute the resulting grade
                unsigned int gradeOfResult = xorIndexToGradeAndPosition[productResultXorIndex(indexMv1,indexMv2)].first;

                if((coefficient !=0) && (gradeOfResult!=gradeOuter) && (gradeOfResult!=gradeInner)){
                    // fill the resulting (gradeMv2-gradeMv1)-vector which is expressed in the orthogonal basis
                    currentProduct.indexOfMv1 =  xorIndexToGradeAndPosition[indexMv1].second;
                    currentProduct.indexOfMv2 =  xorIndexToGradeAndPosition[indexMv2].second;
                    currentProduct.indexOfMv3 = indexOfResult;
                    currentProduct.coefficient= coefficient;
                    // add this product to the fine component of mv3
                    outputListOfProducts[gradeOfResult].push_back(currentProduct);
                }
        }

    return outputListOfProducts;
}


std::vector<std::list<productComponent<double>>> ProductTools::generateExplicitGeometricProductList(const unsigned int gradeMv1,const unsigned int gradeMv2,
                                                                                    const Eigen::SparseMatrix<double, Eigen::ColMajor>& transformationMatrixMv1,
                                                                                    const Eigen::SparseMatrix<double, Eigen::ColMajor>& transformationMatrixMv2,
                                                                                    const std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor> >& transformationMatricesMv3,
                                                                                    const Eigen::VectorXd& diagonalMetric) const{
    std::vector<std::list<productComponent<double>>> outputListOfProducts(dimension+1);
    // As we ignore the inner and outer products contributions. Thus, all products whose resulting grade is gradeOuter=gradeMv1+gradeMv2 or gradeInner=abs(gradeMv1-gradeMv2) will be ignored
    const unsigned int gradeOuter = gradeMv1+gradeMv2;
    const auto gradeInner = (unsigned int)std::abs((float)gradeMv1-gradeMv2);

    // for each element of grade 'gradeMv1'
    for(unsigned int indexMv1=0; indexMv1<(unsigned int)transformationMatrixMv1.rows(); ++indexMv1){
        // for each element of grade 'gradeMv2'
        for(unsigned int indexMv2=0; indexMv2<(unsigned int)transformationMatrixMv2.rows(); ++indexMv2){
            // Might result in higher than 1 resulting grade, thus we create a map of resulting vectors (usage: mv3[grade]= Sparsevector)
            std::map<unsigned int, Eigen::SparseVector<double,Eigen::ColMajor> > mv3Ortho;

            // compute the coefficient of the current product, in the orthogonal basis
            // -> for each basis vector in the orthogonal space of grade 'gradeMv1'
            for(Eigen::SparseMatrix<double>::InnerIterator itMv1(transformationMatrixMv1,indexMv1); itMv1; ++itMv1){
                // -> for each basis vector in the orthogonal space of grade 'gradeMv2'
                for(Eigen::SparseMatrix<double>::InnerIterator itMv2(transformationMatrixMv2,indexMv2); itMv2; ++itMv2){

                    // get the xor indices of Mv1 and Mv2 (independent of the basis)
                    unsigned int xorIndexMv1 = gradePositionToXorIndex[gradeMv1][itMv1.index()];
                    unsigned int xorIndexMv2 = gradePositionToXorIndex[gradeMv2][itMv2.index()];

                    // compute the coefficient
                    double coefficient = productCoefficientFromMetric(xorIndexMv1, xorIndexMv2,diagonalMetric);

                    // compute the resulting index
                    unsigned int indexOfResult = xorIndexToGradeAndPosition[productResultXorIndex(xorIndexMv1, xorIndexMv2)].second;
                    // compute the resulting grade
                    unsigned int gradeOfResult = xorIndexToGradeAndPosition[productResultXorIndex(xorIndexMv1,xorIndexMv2)].first;

                    // fill the resulting (gradeMv2-gradeMv1)-vector which is expressed in the orthogonal basis
                    double resultingCoefficient = coefficient * itMv1.value() * itMv2.value();

                    // in a sparse vector, add only non-zero values, furthermore ignore the outer and inner product contributions
                    if((resultingCoefficient !=0) && (gradeOfResult!=gradeOuter) && (gradeOfResult!=gradeInner)){
                        // If the vector whose grade is gradeOfResult is not already initialized then initialize with the correct size
                        auto itMv3Ortho = mv3Ortho.try_emplace(gradeOfResult, transformationMatricesMv3[gradeOfResult].cols()).first;
                        itMv3Ortho->second.coeffRef(indexOfResult) += resultingCoefficient;
                    }
                }
            } // end of: compute the coefficient of the current product, in the orthogonal basis

            // result in the original basis
            // note that the result might have more than 1 grade thus we have loop over all possible grades of the result
            // except the grades corresponding to the outer and inner products
            // for each grade of the resulting multivector, compute the result in the original basis
            for(const auto & [gradeMv3, vectorMv3Ortho] : mv3Ortho){
                Eigen::SparseVector<double, Eigen::ColMajor> mv3 = transformationMatricesMv3[gradeMv3] * vectorMv3Ortho;
                // add all non-zero component to the resulting list of product component
                for (Eigen::SparseVector<double>::InnerIterator itMv3(mv3); itMv3; ++itMv3) {
                    productComponent<double> currentProduct;
                    // As we ignore the contribution of the outer product and the inner product, we have to check that the grade of the result
                    // is different of (gradeMv1+gradeMv2)  and (gradeMv2-gradeMv1)
                    if (itMv3.value() != 0.0){      // maybe useless test ...
                        currentProduct.indexOfMv1 = indexMv1;
                        currentProduct.indexOfMv2 = indexMv2;
                        currentProduct.indexOfMv3 = itMv3.index();
                        currentProduct.coefficient = itMv3.value();
                        outputListOfProducts[gradeMv3].push_back(currentProduct);
                    }
                }
            }// end of: compute the coefficient of the current product, in the original basis
        }
    }
    return outputListOfProducts;
}


unsigned int ProductTools::getGrade(unsigned int xorIndex) const{
    return xorIndexToGradeAndPosition[xorIndex].first;

}

unsigned int ProductTools::getHomogeneousIndex(unsigned int xorIndex) const{
    return xorIndexToGradeAndPosition[xorIndex].second;
}


unsigned int ProductTools::getXorIndex(unsigned int grade, unsigned int homogeneousIndex) const{
    return gradePositionToXorIndex[grade][homogeneousIndex];
}


double getScaleInversePseudoScalar(const Eigen::MatrixXd &metric){
    // the pseudo scalar I = e1^e2^...^en satisfies reverse(I) I = det(metric), thus
    // I^{-1} = reverse(I) / det(metric) = (-1)^{n(n-1)/2} I / det(metric).
    // This function returns 1/det(metric), the caller handles the reverse sign.
    // (0 for a degenerate metric, the pseudo scalar is then not invertible)
    const double determinant = metric.determinant();
    if(determinant == 0.0)
        return 0.0;
    return 1.0 / determinant;
}
