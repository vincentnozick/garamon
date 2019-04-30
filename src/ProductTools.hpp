// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// ProductTools.hpp
// This file is part of the Garamon Generator.
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

/// \file ProductTools.hpp
/// \author Vincent Nozick, Stephane Breuils 
/// \brief Geometric algebra "naive" implementation tools dedicated to build the optimized C++ products.


#ifndef GARAGEN_PRODUCTTOOLS_HPP
#define GARAGEN_PRODUCTTOOLS_HPP

#include <vector>
#include <array>
#include <map>
#include <utility>     // pairs
#include <Eigen/Dense> // use the metric which is dense
#include <Eigen/Sparse> // use the transformation matrices which are sparse
#include <iostream>
#if defined(_MSC_BUILD)
	#ifndef __builtin_popcount
		#define __builtin_popcount __popcnt
	#endif
#endif // __WINDOWS MSVC compiler__


/// In the computation of mv3 = mv1^mv2, this represents a quadruplet containing the indices of mv1,mv2,mv3 and a coefficient in a product between two blades.
/// e.g. suppose that a product is defined by mv3[2] += -3.0*mv1[3]*mv2[4]  then the first quadruplet associated with this product will be 3,4,1,-3.0
template<typename T>
struct productComponent{
    unsigned int indexOfMv1; // in the example above, indexOfMv1 would be 3
    unsigned int indexOfMv2; // in the example above, indexOfMv2 would be 4
    unsigned int indexOfMv3; // in the example above, indexOfMv2 would be 4
    T coefficient;           // in the example above, coefficient would be -1
};

template<typename C>
bool compareProductComponents(const C &a, const C &b)
//bool operator < (const productComponent &a, const productComponent &b) const
{
    // first compare on indexOfMv3
    if(a.indexOfMv3 > b.indexOfMv3) return false;
    if(a.indexOfMv3 < b.indexOfMv3) return true;

    // second compare on indexOfMv1
    if(a.indexOfMv1 > b.indexOfMv1) return false;
    if(a.indexOfMv1 < b.indexOfMv1) return true;

    // thrid compare on indexOfMv2
    if(a.indexOfMv2 > b.indexOfMv2) return false;
    if(a.indexOfMv2 < b.indexOfMv2) return true;

    return true;
}


class ProductTools {
public:

    /// constructor
    /// \param vectorSpaceDimension is the dimension of the vector space of the algebra
    ProductTools(const unsigned int vectorSpaceDimension);

    /// destructor
    ~ProductTools();

    /// \brief compute the Hamming weight of the xor index (blade) xorIndexMv (i.e the number of bits=1 in the binary value xorIndexMv)
    /// \param xorIndexMv is a xorIndex, i.e. a each of his bit refers to a basis blade of the vector space of the algebra
    static unsigned int hammingWeight(const unsigned int xorIndexMv);

    /// \brief return whether the outer product between the two blades whose indices are xorIndexMv1 and xorIndexMv2 is a null blade
    /// \param xorIndexMv1 is a xorIndex, i.e. a each of his bit refers to a basis blade of the vector space of the algebra
    /// \param xorIndexMv2 is a xorIndex, i.e. a each of his bit refers to a basis blade of the vector space of the algebra
    static bool outerProductExists(const unsigned int xorIndexMv1, const unsigned int xorIndexMv2);

    /// \brief return the sign of the outer product between the two blades whose indices are xorIndexMv1 and xorIndexMv2
    /// \param xorIndexMv1 is a xorIndex, i.e. a each of his bit refers to a basis blade of the vector space of the algebra
    /// \param xorIndexMv2 is a xorIndex, i.e. a each of his bit refers to a basis blade of the vector space of the algebra
    static int outerProductSign(const unsigned int xorIndexMv1, const unsigned int xorIndexMv2);

    /// \brief return whether the inner product exists
    /// \param xorIndexMv1 is a xorIndex, i.e. a each of his bit refers to a basis blade of the vector space of the algebra
    /// \param xorIndexMv2 is a xorIndex, i.e. a each of his bit refers to a basis blade of the vector space of the algebra
    static bool innerProductExists(const unsigned int xorIndexMv1, const unsigned int xorIndexMv2);


    /// \brief Compute the coefficient required in the computation of the inner and geometric products of 2 multivectors.
    /// \param xorIndexMv1 is a xorIndex, i.e. a each of his bit refers to a basis blade of the vector space of the algebr  in the algebra related orthogonal space
    /// \param xorIndexMv2 is a xorIndex, i.e. a each of his bit refers to a basis blade of the vector space of the algebra in the algebra related orthogonal space
    /// \param diagonalMetric defines the metric in the algebra related orthogonal space
    static double productCoefficientFromMetric(const unsigned int xorIndexMv1, const unsigned int xorIndexMv2,
                                               const Eigen::VectorXd &diagonalMetric);

    /// \brief Compute the position of the product between two XOR indices. Is common for the geometric, outer and inner product
    /// \param xorIndexMv1 is a xorIndex, i.e. a each of his bit refers to a basis blade of the vector space of the algebra
    /// \param xorIndexMv2 is a xorIndex, i.e. a each of his bit refers to a basis blade of the vector space of the algebra
    static unsigned int productResultXorIndex(const unsigned int xorIndexMv1, const unsigned int xorIndexMv2);

    /// \brief Generate all the components of the product mv1{gradeMv1} ^ mv2{gradeMv2}.
    /// \param xorIndexMv1 is a xorIndex, i.e. a each of his bit refers to a basis blade of the vector space of the algebr  in the algebra related orthogonal space
    /// \param xorIndexMv2 is a xorIndex, i.e. a each of his bit refers to a basis blade of the vector space of the algebra in the algebra related orthogonal space
    /// \return For each component of mv3=mv1^mv2 with grade gradeMv3=gradeMv1+gradeMv2, we return a list of quadruplet (indexMv1,indexMv2,indexMv3,coefficient).
    std::list<productComponent<double>> generateExplicitOuterProductList(const unsigned int gradeMv1, const unsigned int gradeMv2) const;


    /// \brief Generate all the components of the inner product mv1<gradeMv1> . mv2<gradeMv2> in a Euclidean space with diagonal metric. In practice, we compute here the Hestenes inner product.
    /// \param gradeMv1 the grade of the homogeneous multivector Mv1 to be considered
    /// \param gradeMv2 the grade of the homogeneous multivector Mv2 to be considered
    /// \param diagonalMetric defines the metric in the algebra related orthogonal space
    /// \return for each component of the result mv3=mv1.mv2 with grade gradeMv3=|gradeMv1-gradeMv2|, we return a list of quadruplet (indexMv1,indexMv2,indexMv3,coefficient).
    std::list<productComponent<double>> generateExplicitInnerProductListEuclideanSpace(const unsigned int gradeMv1, 
                                                                                                     const unsigned int gradeMv2,
                                                                                                     const Eigen::VectorXd &diagonalMetric) const;

    /// \brief Generate all the components of the inner product mv1<gradeMv1> . mv2<gradeMv2>. In practice, we compute here the Hestenes inner product.
    /// \param gradeMv1 the grade of the homogeneous multivector Mv1 to be considered
    /// \param gradeMv2 the grade of the homogeneous multivector Mv2 to be considered
    /// \param transformationMatrixMv1 is the matrix that put a k-vector of grade "gradeMv1" to the algebra related orthogonal space
    /// \param transformationMatrixMv2 is the matrix that put a k-vector of grade "gradeMv2" to the algebra related orthogonal space
    /// \param transformationMatrixMv3 is the matrix that put back a k-vector of grade "gradeMv3" from the algebra related orthogonal space to the original space
    /// \param diagonalMetric defines the metric in the algebra related orthogonal space
    /// \return for each component of the result mv3=mv1.mv2 with grade gradeMv3=|gradeMv1-gradeMv2|, we return a list of quadruplet (indexMv1,indexMv2,indexMv3,coefficient).
    std::list<productComponent<double>> generateExplicitInnerProductList(const unsigned int gradeMv1,
                                                                         const unsigned int gradeMv2,
                                                                         const Eigen::SparseMatrix<double, 0> &transformationMatrixMv1,
                                                                         const Eigen::SparseMatrix<double, 0> &transformationMatrixMv2,
                                                                         const Eigen::SparseMatrix<double, 0> &transformationMatrixMv3,
                                                                         const Eigen::VectorXd &diagonalMetric) const;



    /// \brief Generate all the components of the geometric product mv1<gradeMv1> mv2<gradeMv2>.
    /// \param gradeMv1 the grade of the homogeneous multivector Mv1 to be considered
    /// \param gradeMv2 the grade of the homogeneous multivector Mv2 to be considered
    /// \param transformationMatrixMv1 is the matrix that put a k-vector of grade "gradeMv1" to the algebra related orthogonal space
    /// \param transformationMatrixMv2 is the matrix that put a k-vector of grade "gradeMv2" to the algebra related orthogonal space
    /// \param transformationMatrixMv3 is the matrix that put back a k-vector of grade "gradeMv3" from the algebra related orthogonal space to the original space
    /// \param diagonalMetric defines the metric in the algebra related orthogonal space
    /// \return a vector of products (each cell of the vector contains its associated grade products). For each component of the result mv3=mv1 mv2, we return a list of quadruplet (indexMv1,indexMv2,indexMv3,coefficient).
    std::vector<std::list<productComponent<double>>> generateExplicitGeometricProductListEuclideanSpace(const unsigned int gradeMv1,const unsigned int gradeMv2,
                                                                                          const Eigen::VectorXd& diagonalMetric) const;



    /// \brief Generate all the components of the geometric product mv1<gradeMv1> mv2<gradeMv2>.
    /// \param gradeMv1 the grade of the homogeneous multivector Mv1 to be considered
    /// \param gradeMv2 the grade of the homogeneous multivector Mv2 to be considered
    /// \param transformationMatrixMv1 is the matrix that put a k-vector of grade "gradeMv1" to the algebra related orthogonal space
    /// \param transformationMatrixMv2 is the matrix that put a k-vector of grade "gradeMv2" to the algebra related orthogonal space
    /// \param transformationMatrixMv3 is the matrix that put back a k-vector of grade "gradeMv3" from the algebra related orthogonal space to the original space
    /// \param diagonalMetric defines the metric in the algebra related orthogonal space
    /// \return a vector of products (each cell of the vector contains its associated grade products). For each component of the result mv3=mv1 mv2, we return a list of quadruplet (indexMv1,indexMv2,indexMv3,coefficient).
    std::vector<std::list<productComponent<double>>> generateExplicitGeometricProductList(const unsigned int gradeMv1,const unsigned int gradeMv2,
                                                                                          Eigen::SparseMatrix<double, Eigen::ColMajor>& transformationMatrixMv1,
                                                                                          Eigen::SparseMatrix<double, Eigen::ColMajor>& transformationMatrixMv2,
                                                                                          std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor> >& transformationMatricesMv3,
                                                                                          const Eigen::VectorXd& diagonalMetric) const;



    // tool to debug
    void displayGradeToXorIndices() const;

    /// \brief getter: return the grade knowing the Xor index
    unsigned int getGrade(unsigned int xorIndex) const;

    /// \brief getter: return the homogeneous index knowing the Xor index
    unsigned int getHomogeneousIndex(unsigned int xorIndex) const;

    /// \brief getter: return the xor index knowing the grade and the homogeneous index
    unsigned int getXorIndex(unsigned int grade, unsigned int homogeneousIndex) const;


protected:
    const unsigned int dimension;

    /// order induced by the xor:                 scal, 1, 2, 12, 3, 13, 23, 123
    /// Mapping between
    /// Xor index:                                   0, 1, 2,  3, 4,  5,  6, 7
    /// and grade and position in the sequence:      0, 0, 1,  0, 2,  1,  2, 0

    // defines the one to one correspondence between a xor index (index in the multivector) and the grade and position in the homogeneous vector
    // grade and position to xor index
    /// xorIndex = gradePositionToXorIndex[grade][pos];
    std::vector<std::vector<unsigned int> > gradePositionToXorIndex;

    /// array that converts xor index to grade and position
    /// first: grade
    /// second: position
    std::vector<std::pair<unsigned int,unsigned int> > xorIndexToGradeAndPosition; // number of pairs is 2^dimension
};

/// \brief generate a scale containing the sign of the norm of the pseudo-inverse
double getScaleInversePseudoScalar(const Eigen::MatrixXd &metric);

#endif //GARAGEN_PRODUCTTOOLS_HPP
