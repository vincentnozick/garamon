// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// Mvec.hpp
// This file is part of the Garamon for project_namespace.
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

/// \file Mvec.hpp
/// \author Stephane Breuils, Vincent Nozick
/// \brief Class to define a multivector and its basic operators in the Geometric algebra of project_namespace.


// Anti-doublon
#ifndef project_inclusion_guard
#define project_inclusion_guard
#pragma once

// External Includes
#include <Eigen/Core>
#include <algorithm>
#include <list>
#include <iostream>
#include <cmath>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

// Internal Includes (Mvec.hpp is the header to include to use the library: it exports the constants and the utilities)
#include "project_namespace/Utility.hpp"   // IWYU pragma: export
#include "project_namespace/Constants.hpp" // IWYU pragma: export

#include "project_namespace/Outer.hpp"
#include "project_namespace/Inner.hpp"
#include "project_namespace/Geometric.hpp"

#include "project_namespace/OuterExplicit.hpp"
#include "project_namespace/InnerExplicit.hpp"
#include "project_namespace/GeometricExplicit.hpp"

/*!
 * @namespace project_namespace
 */
namespace project_namespace{

    /// \brief the scalar types a multivector can be constructed from or combined with (e.g. 2.0 * mv):
    /// the arithmetic types only, so that the implicit constructor Mvec(S) does not make any other type
    /// (e.g. a char *) convertible to a multivector
    template<typename S>
    concept Scalar = std::is_arithmetic_v<S>;


    /// \class Kvec
    /// \brief class defining a single grade component of a multivector.
    template<class T>
    struct Kvec{

        Eigen::Matrix<T, Eigen::Dynamic, 1> vec;  /*!< dynamic vector of Eigen Library */

        unsigned int grade; /*!< grade k of the k-vector */

        /// \brief operator == to test the equality between 2 k-vectors.
        bool operator == (const Kvec<T>& other) const {
            return vec == other.vec;
        }
    };



    /// \class Mvec
    /// \brief class defining multivectors.
    template<typename T = double>
    class Mvec {

    protected:
        std::list<Kvec<T>> mvData;  /*!< set of k-vectors, mapped by grade */
        unsigned int gradeBitmap;   /*!< ith bit to 1 if grade i is contained in the multivector */

    public:

        /// \brief Default constructor, generate an empty multivector equivalent to the scalar 0.
        Mvec();

        /// \brief Copy constructor
        /// \param mv - the multivector to be copied
        Mvec(const Mvec& mv);

        /// \brief Copy constructor of Mvec
        /// \param mv - the multivector which has to be copied
        Mvec(Mvec&& mv); // move constructor


        template <typename U>
        friend class Mvec;

        /// \brief Copy constructor of Mvec from different types
        /// \param mv - the multivector with another template type
        template<typename U>
        Mvec(const Mvec<U> &mv);

        /// \brief Constructor of Mvec from a scalar
        /// \param val - scalar value
        template<Scalar S>
        Mvec(const S val);

        /// Destructor
        ~Mvec();

        /// \brief Overload the assignment operator. No need any copy when the argument is an R-value, just need to move.
        /// \param mv - Mvec used as a Rvalue
        /// \return assign other to this
        Mvec& operator=(Mvec&& mv);

        /// \brief Overload the assignment operator
        /// \param mv - Mvec
        /// \return assign mv to this object
        Mvec& operator=(const Mvec& mv);

        /// \brief defines the addition between two Mvec
        /// \param mv2 - second operand of type Mvec
        /// \return this + mv2
        Mvec operator+(const Mvec &mv2) const;

        /// \brief defines the addition between a Mvec and a scalar
        /// \param value - second operand (scalar)
        /// \return this + scalar
        template<Scalar S>
        Mvec operator+(const S &value) const;

        /// \brief defines the addition between a scalar and a Mvec
        /// \param value a scalar
        /// \param mv - second operand of type Mvec
        /// \return scalar + mv
        template<typename U, Scalar S>
        friend Mvec<U> operator+(const S &value, const Mvec<U> &mv);

        /// \brief Overload the += operator, corresponds to this += mv
        /// \param mv - Mvec to be added to this object
        /// \return this += mv
        Mvec& operator+=(const Mvec& mv);

        /// \brief defines the opposite of a multivector
        /// \param mv: operand of type Mvec
        /// \return -mv
        template<typename U>
        friend Mvec<U> operator-(const Mvec<U> &mv); // unary operator -mv1

        /// \brief defines the difference between two Mvec
        /// \param mv2 - second operand of type Mvec
        /// \return this - mv2
        Mvec<T> operator-(const Mvec<T> &mv2) const;

        /// \brief defines the difference between a Mvec and a scalar
        /// \param value - second operand (scalar)
        /// \return this - scalar
        template<Scalar S>
        Mvec operator-(const S &value) const;

        /// \brief defines the difference between a scalar and a Mvec
        /// \param value a scalar
        /// \param mv - second operand of type Mvec
        /// \return scalar - mvX
        template<typename U, Scalar S>
        friend Mvec<U> operator-(const S &value, const Mvec<U> &mv);

        /// \brief Overload the -= operator, corresponds to this -= mv
        /// \param mv - Mvec to be added to this object
        /// \return this -= mv
        Mvec& operator-=(const Mvec& mv);

        /// \brief defines the outer product between two multivectors
        /// \param mv2 - a multivector
        /// \return this^mv2
        Mvec operator^(const Mvec &mv2) const;

        /// \brief defines the outer product a multivector and a scalar
        /// \param value - a scalar
        /// \return this^value
        template<Scalar S>
        Mvec operator^(const S &value) const;

        /// \brief defines the outer product between a scalar and a multivector
        /// \param value - a scalar
        /// \param mv - a multivector
        /// \return value ^ mv
        template<typename U, Scalar S>
        friend Mvec<U> operator^(const S &value, const Mvec<U> &mv);

        /// \brief Overload the outer product with operator =, corresponds to this ^= mv
        /// \param mv - Mvec to be wedged to this object
        /// \return this ^= mv
        Mvec& operator^=(const Mvec& mv);

        /// \brief defines the inner product between two multivectors
        /// \param mv2 - a multivector
        /// \return this.mv2
        Mvec operator|(const Mvec &mv2) const;

        /// \brief defines the inner product between a multivector and a scalar
        /// \param value - a scalar
        /// \return this.value = 0 (empty multivector)
        template<Scalar S>
        Mvec operator|(const S &value) const;

        /// \brief defines the inner product between a scalar and a multivector
        /// \param value - a scalar
        /// \param mv - a multivector
        /// \return mv.value = 0 (empty multivector)
        template<typename U, Scalar S>
        friend Mvec<U> operator|(const S &value, const Mvec<U> &mv);

        /// \brief Overload the inner product with operator =, corresponds to this |= mv
        /// \param mv - Mvec to be computed in the inner product with this object
        /// \return this |= mv
        Mvec& operator|=(const Mvec& mv);

        /// \brief defines the right contraction between two multivectors
        /// \param mv2 - a multivector
        /// \return the right contraction : $this \\lfloor mv2$
        Mvec operator>(const Mvec &mv2) const;

        /// \brief defines the right contraction between a multivector and a scalar
        /// \param value - a scalar
        /// \return $this \\lfloor value$
        template<Scalar S>
        Mvec operator>(const S &value) const;

        /// \brief defines the right contraction between a scalar and a multivector
        /// \param value - a scalar
        /// \param mv - a multivector
        /// \return $value \\lfloor this = 0$ (empty multivector)
        template<typename U, Scalar S>
        friend Mvec<U> operator>(const S &value, const Mvec<U> &mv);

        /// \brief defines the left contraction between two multivectors
        /// \param mv2 - a multivector
        /// \return the left contraction $this \\rfloor mv2$
        Mvec operator<(const Mvec &mv2) const;

        /// \brief defines the left contraction between a multivector and a scalar
        /// \param value - a scalar
        /// \return $value \\rfloor this$ = 0 (empty multivector)
        template<Scalar S>
        Mvec operator<(const S &value) const;

        /// \brief defines the left contraction between a scalar and a multivector
        /// \param value - a scalar
        /// \param mv - a multivector
        /// \return $this \\rfloor value$
        template<typename U, Scalar S>
        friend Mvec<U> operator<(const S &value, const Mvec<U> &mv);

        /// \brief defines the geometric product between two multivectors
        /// \param mv2 - a multivector
        /// \return mv1*mv2
        Mvec operator*(const Mvec &mv2) const;

    project_singular_metric_comment_begin
        /// \brief defines the outer product between two multivectors, where the second multivector is dualized during the product computation.
        /// \param mv2 - a multivector that will be dualized during the wedge with the calling multivector.
        /// \return a multivector.
        Mvec<T> outerPrimalDual(const Mvec<T> &mv2) const;

        /// \brief defines the outer product between two multivectors, where the first multivector is dualized during the product computation.
        /// \param mv2 - a primal form of a multivector; the object multivector will be dualized during the wedge with the calling multivector.
        /// \return a multivector.
        Mvec<T> outerDualPrimal(const Mvec<T> &mv2) const;

        /// \brief defines the outer product between two multivectors, where both multivectors are dualized during the product computation.
        /// \param mv2 - a primal form of a multivector; the object multivector will be dualized during the wedge with the calling multivector.
        /// \return a multivector.
        Mvec<T> outerDualDual(const Mvec<T> &mv2) const;
    project_singular_metric_comment_end
        

        /// \brief defines the scalar product between two multivectors (sum of the inner products between same grade pairs from the 2 multivectors)
        /// \param mv2 - a multivector
        /// \return a scalar
        Mvec<T> scalarProduct(const Mvec<T> &mv2) const;

        /// \brief defines the dot product between two multivectors (Inner product + scalar product)
        /// \param mv2 - a multivector
        /// \return a multivector.
        Mvec<T> dotProduct(const Mvec<T> &mv2) const;

        /// \brief defines the geometric product between a multivector and a scalar
        /// \param value - a scalar
        /// \return mv2*value
        template<Scalar S>
        Mvec operator*(const S &value) const;

        /// \brief defines the geometric product between a scalar and a multivector
        /// \param value - a scalar
        /// \param mv - a multivector
        /// \return value*mv2
        template<typename U, Scalar S>
        friend Mvec<U> operator*(const S &value, const Mvec<U> &mv);

        /// \brief Overload the geometric product with operator =, corresponds to this *= mv
        /// \param mv - Mvec to be multiplied to this object
        /// \return this *= mv
        Mvec& operator*=(const Mvec& mv);

        /// \brief defines the geometric product with a multivector and the inverse of a second multivector
        /// \param mv2 - a multivector
        /// \return this / mv2
        Mvec operator/(const Mvec &mv2) const;

        /// \brief defines the scalar inverse between a multivector and a scalar
        /// \param value - a scalar
        /// \return this / value
        template<Scalar S>
        Mvec operator/(const S &value) const;

        /// \brief defines the inverse product between a scalar and a multivector
        /// \param value - a scalar
        /// \param mv - a multivector
        /// \return value / mv
        template<typename U, Scalar S>
        friend Mvec<U> operator/(const S &value, const Mvec<U> &mv);

        /// \brief Overload the inverse with operator =, corresponds to this /= mv
        /// \param mv - Mvec to be inverted to this object
        /// \return this /= mv
        Mvec& operator/=(const Mvec& mv);


        /// \brief Overload the inverse with operator =, corresponds to this /= value
        /// \param value - scalar to divide the object
        /// \return this /= value
        template<Scalar S>
        Mvec& operator/=(const S &value);


        /// \brief the reverse of a multivector, i.e. if mv = a1^a2^...^an, then reverse(mv) = an^...^a2^a1
        /// \param mv - a multivector
        /// \return reverse of mv
        template<typename U>
        friend Mvec<U> operator~(const Mvec<U> &mv);

        /// \brief the dual of a k-vector is defined as $A_k^* = A_k \\lcont I_n$ (left contraction onto the pseudo-scalar), for a multivector, we just dualize all its components. If the metric is degenerated, this function computes the right complement (mv ^ !mv = I).
        /// \param mv - a multivector
        /// \return the dual of mv
        template<typename U>
        friend Mvec<U> operator!(const Mvec<U> &mv);

        /// \brief boolean operator that tests the equality between two Mvec (mv1 != mv2 is rewritten by the compiler as !(mv1 == mv2))
        /// \param mv2 - second operand of type Mvec
        /// \return whether two Mvec have the same coefficients
        inline bool operator==(const Mvec& mv2) const {
            if(gradeBitmap != mv2.gradeBitmap)
                return false;

            return mvData == mv2.mvData;   //// list : ca marche que si les listes sont ordonnees
        }

        /// \brief compute the inverse of a multivector
        /// \return - the inverse of the current multivector
        Mvec<T> inv() const;

        /// \brief Display all the non-null basis blades of this objects
        /// \param stream - destination stream
        /// \param mvec - the multivector to be outputed
        /// \return a stream that contains the list of the non-zero element of the multivector
        template<typename U>
        friend std::ostream& operator<< (std::ostream& stream, const Mvec<U> &mvec);

        /// \brief overload the casting operator, using this function, it is now possible to compute : float a = float(mv);
        /// (explicit: a multivector is not silently converted to its scalar part, e.g. in fabs(mv) or double d = mv;
        /// a template, so that a cast to any arithmetic type, e.g. int(mv), works as before)
        /// \return the scalar part of the multivector
        template<Scalar S>
        explicit operator S () const {
            // if no scalar part, return

            if( (gradeBitmap & 1) == 0 )
                return S(0);

            // assuming now that the scalar part exists, return it
            return S(findGrade(0)->vec[0]);
        }

        /// \brief Overload the [] operator to assign a basis blade to a multivector. As an example, float a = mv[E12] = 42.
        /// \param idx - the basis vector index related to the query.
        /// \return the coefficient of the multivector corresponding to the "idx" component.
        /// \todo handle the cases when the number of indices is higher than the dimension, and when the index is too high
        T& operator[](const int idx ){
            // the k-vector is created (null) if the multivector does not contain it yet
            return createVectorXdIfDoesNotExist(xorIndexToGrade[idx])->vec[xorIndexToHomogeneousIndex[idx]];
        }

        /// \brief Overload the [] operator to copy a basis blade of this multivector. As an example, float a = mv[E12].
        /// \param idx - the basis vector index related to the query.
        /// \return the coefficient of the multivector corresponding to the "idx" component.
        const T& operator[](const int idx) const{
            // if the multivector does not contain the k-vector, return a reference on zero
            auto it = findGrade(xorIndexToGrade[idx]);
            return (it == mvData.end()) ? zero<T> : it->vec[xorIndexToHomogeneousIndex[idx]];
        }

        /// \cond DEV
        /// \brief returns a multivector with one component to 1
        /// \param grade : grade of the component to enable
        /// \param index : index of the parameter in the k-vector (k = grade)
        inline Mvec componentToOne(const unsigned int grade, const int index){
            Mvec mv1;
            mv1.mvData.push_back(Kvec<T>{.vec = Eigen::Matrix<T, Eigen::Dynamic, 1>::Unit(binomialArray[grade], index), .grade = grade});
            mv1.gradeBitmap = 1u << grade;

            return mv1;
        }
        /// \endcond // do not comment this functions

        /// \cond DEV
        /// \brief create a VectorXd if it has not yet been created
        /// \param grade - grade of the considered kvector
        /// \return nothing
        inline typename std::list<Kvec<T>>::iterator createVectorXdIfDoesNotExist(const unsigned int grade){

            // the k-vectors are sorted by increasing grade: find the first one whose grade is not lower than 'grade'
            auto it = std::ranges::find_if(mvData, [grade](const unsigned int kvecGrade){ return kvecGrade >= grade; }, &Kvec<T>::grade);
            if(it != mvData.end() && it->grade == grade)
                return it;

            // not found: insert a null k-vector before the current element (possibly the end)
            gradeBitmap |= 1u << grade;
            return mvData.insert(it, Kvec<T>{.vec = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(binomialArray[grade]), .grade = grade});
        }


        /// \cond DEV
        /// \brief remove the k-vector pointed by 'it' if all its components are 0: the result of an operation never contains null k-vectors (the empty multivector is 0).
        /// \param it - iterator on a k-vector of this multivector
        inline void eraseIfNull(typename std::list<Kvec<T>>::iterator it){
            if(!((it->vec.array() != T(0)).any())){
                gradeBitmap &= ~(1u << it->grade);
                mvData.erase(it);
            }
        }
        /// \endcond // do not comment this functions

        /// \cond DEV
        /// \brief modify the element of the multivector whose grade is "grade"
        /// \param grade : the grade to access
        /// \param indexVectorXd : index of the k-vector (with k = "grade")
        /// \return the element of the Mv whose grade is grade and index in the VectorXd is indexVectorXd
        inline T& at(const int grade, const int indexVectorXd){
            auto it = createVectorXdIfDoesNotExist(grade);
            return it->vec.coeffRef(indexVectorXd);  
        }
        /// \endcond // do not comment this functions

        /// \cond DEV
        /// \brief return the element of the multivector whose grade is "grade"
        /// \param grade : the grade to access
        /// \param indexVectorXd : index of the k-vector (with k = "grade")
        /// \return the element of the Mv whose grade is grade and index in the VectorXd is indexVectorXd
        inline T at(const int grade, const int indexVectorXd) const{

            if((gradeBitmap & (1<<(grade))) == 0 )
                return T(0);

            auto it = findGrade(grade);
            return it->vec.coeff(indexVectorXd);
        }
        /// \endcond // do not comment this functions

        /// \brief the L2-norm of the mv is sqrt( abs( mv.mv ) )
        /// \return the L2-norm of the multivector (as a double)
        T inline norm() const {
            return sqrt( fabs( T((*this).scalarProduct( this->reverse() )) ));
        }

        /// \brief the L2-norm over 2 of the mv is mv.mv
        /// \return the L2-norm of the multivector (as a double)
        T inline quadraticNorm() const {
            return T( this->reverse() > (*this) );
        }

        /// \brief compute the dual of a multivector (i.e. mv* = mv < I, the left contraction of mv onto the pseudo-scalar). If the metric is degenerated, this function computes the right complement (mv ^ !mv = I).
        /// \return - the dual of the multivector
        Mvec<T> dual() const;

        /// \brief compute the reverse of a multivector
        /// \return - the reverse of the multivector
        Mvec<T> reverse() const;

        /// \brief search in the multivector for a Kvec of grade "grade"
        /// \return return a const iterator on the Kvec if exist, else return mvData.end()
        inline typename std::list<Kvec<T>>::const_iterator findGrade(const unsigned int & gradeToFind) const {
            return std::ranges::find(mvData, gradeToFind, &Kvec<T>::grade);
        }

        /// \brief search in the multivector for a Kvec of grade "grade"
        /// \return return an iterator on the Kvec if exist, else return mvData.end()
        inline typename std::list<Kvec<T>>::iterator findGrade(const unsigned int & gradeToFind) {
            return std::ranges::find(mvData, gradeToFind, &Kvec<T>::grade);
        }

        /// \brief return the (highest) grade of the multivector
        /// \return the highest grade of the multivector
        inline int grade() const {
            return (mvData.begin() == mvData.end()) ? 0 : mvData.rbegin()->grade;
        }

        /// \brief return the all non-zero grades of the multivector (several grades for non-homogeneous multivectors)
        /// \return a list of the grades encountered in the multivector
        std::vector<unsigned int> grades() const;

        /// \brief returns a multivector that contains all the components of this multivector whose grade is i
        /// \return an empty Mvec if the requested element is not part of the multivector, or the multivector that contains only this element if present in the current multivector.
        Mvec grade(const int i) const;

        /// \brief tell whether the multivector has grade component
        /// \param grade - grade of the considered kvector
        /// \return true if multivector has grade component, false else
        inline bool isGrade(const unsigned int grade) const{
            if((grade == 0) && isEmpty()) return true;
            return ( (gradeBitmap & (1<<grade)) != 0 );
        }

        /// \brief partially or completely erase the content of a multivector
        /// \param grade : if < 0, erase all the multivector, else just erase the part of grade "grade".
        void clear(const int grade = -1);

        /// \brief check is a mutivector is empty, i.e. corresponds to 0.
        /// \return True if the multivector is empty, else False.
        inline bool isEmpty() const {
            return mvData.empty();
        }

        /// \brief A multivector is homogeneous if all its components have the same grade.
        /// \return True if the multivector is homogeneous, else False.
        inline bool isHomogeneous() const {
            return mvData.size() < 2; // only one element, or zero element on the list
        }

        /// \brief inplace simplify the multivector such that all the values with a magnitude lower than a epsilon in the Mv are set to 0.
        /// \param epsilon - threshold, with default value the epsilon of the float/double/long double type from numeric_limits.
        void roundZero(const T epsilon = std::numeric_limits<T>::epsilon());

        /// \brief Specify if two multivectors have the same grade.
        /// \param mv - multivector to compare with.
        /// \return true if the two multivectors have the same grade, else return false.
        inline bool sameGrade(const Mvec<T>& mv) const{
            return grade() == mv.grade();
        }

        /// \brief Display the multivector data (per grade value)
        void display() const;

        /// \cond DEV
        /// \brief a function to output all the element of a k-vector in the multivector indexing order.
        /// \tparam T - type of the multivector
        /// \param stream - stream that will contain the result
        /// \param mvec - multivector to be processed
        /// \param gradeMV - the considered grade
        /// \param moreThanOne - true if it the first element to display (should we put a '+' before)
        template<typename U>
        friend void traverseKVector(std::ostream &stream, const Eigen::Matrix<U, Eigen::Dynamic, 1>, unsigned int gradeMV, bool& moreThanOne);
        /// \endcond // do not comment this functions


        /// \brief functions that enables to extract one component of this multivector and initialize another mv with this component
        /// \tparam Arguments - denotes a variadic list of index
        /// \param grade : the considered grade
        /// \param sizeOfKVector: C(dimension, grade)
        /// \param indexInKvector: the considerd index
        /// \return a new mv with the right component at the right place
        Mvec extractOneComponent(const int grade, const int sizeOfKVector, const int indexInKvector) const {
            Mvec mv;
            auto itThisMV = this->findGrade(grade);
            // null component: the result is the empty multivector
            if(itThisMV==this->mvData.end() || itThisMV->vec.coeff(indexInKvector) == T(0)){
                return mv;
            }
            mv = mv.componentToOne(grade,indexInKvector);
            mv.mvData.begin()->vec.coeffRef(indexInKvector) = itThisMV->vec.coeff(indexInKvector);
            return mv;
        }


        /// \brief set of functions that return multivectors containing only the value '1' for the specified component.
project_multivector_one_component

    };  // end of the class definition


    /* ------------------------------------------------------------------------------------------------ */


    template<typename T>
    Mvec<T>::Mvec():gradeBitmap(0)
    {}


    template<typename T>
    Mvec<T>::Mvec(const Mvec& mv) : mvData(mv.mvData), gradeBitmap(mv.gradeBitmap)
    {
    }


    template<typename T>
    Mvec<T>::Mvec(Mvec<T>&& multivector) : mvData(std::move(multivector.mvData)), gradeBitmap(multivector.gradeBitmap)
    {
    }


    template<typename T>
    template<typename U>
    Mvec<T>::Mvec(const Mvec<U> &mv) : gradeBitmap(mv.gradeBitmap)
    {
        for(const auto & kvec : mv.mvData)
            mvData.push_back(Kvec<T>{.vec = kvec.vec.template cast<T>(), .grade = kvec.grade});
    }


    template<typename T>
    template<Scalar U>
    Mvec<T>::Mvec(const U val) : gradeBitmap(0) {
        if(val != U(0)) {
            gradeBitmap = 1;
            mvData.push_back(Kvec<T>{.vec = Eigen::Matrix<T, Eigen::Dynamic, 1>::Constant(1, T(val)), .grade = 0});
        }
    }


    template<typename T>
    Mvec<T>::~Mvec()
    {}


    template<typename T>
    Mvec<T>& Mvec<T>::operator=(const Mvec& mv){
        if(&mv == this) return *this;
        gradeBitmap = mv.gradeBitmap;
        mvData = mv.mvData;
        return *this;
    }


    template<typename T>
    Mvec<T>& Mvec<T>::operator=(Mvec&& mv){
        mvData = std::move(mv.mvData);
        gradeBitmap = mv.gradeBitmap;
        return *this;
    }


    template<typename T>
    Mvec<T> Mvec<T>::operator+(const Mvec<T> &mv2) const {
        Mvec<T> mv3(*this);
        for(auto & itMv : mv2.mvData) {
            auto it = mv3.createVectorXdIfDoesNotExist(itMv.grade);
            it->vec += itMv.vec;
            mv3.eraseIfNull(it);
        }
        return mv3;
    }


    template<typename U, Scalar S>
    Mvec<U> operator+(const S &value, const Mvec<U> &mv){
        return mv + value;
    }


    template<typename T>
    template<Scalar S>
    Mvec<T> Mvec<T>::operator+(const S &value) const {
        Mvec<T> mv(*this);
        if(value != T(0)) {
            auto it = mv.createVectorXdIfDoesNotExist(0);
            it->vec.coeffRef(0) += value;
            mv.eraseIfNull(it);
        }
        return mv;
    }


    template<typename T>
    Mvec<T>& Mvec<T>::operator+=(const Mvec& mv){
        *this = *this + mv;
        return *this;
    }


    template<typename T>
    Mvec<T> operator-(const Mvec<T> &mv) { // unary -
        Mvec<T> mv2(mv);
        for(auto & itMv : mv2.mvData)
            itMv.vec = -itMv.vec;
        return mv2;
    }


    template<typename T>
    Mvec<T> Mvec<T>::operator-(const Mvec<T> &mv2) const {
        Mvec<T> mv3(*this);
        for(auto & itMv : mv2.mvData) {
            auto it = mv3.createVectorXdIfDoesNotExist(itMv.grade);
            it->vec -= itMv.vec;
            mv3.eraseIfNull(it);
        }
        return mv3;
    }


    template<typename T>
    Mvec<T>& Mvec<T>::operator-=(const Mvec& mv){
        *this = *this - mv;
        return *this;
    }


    template<typename U, Scalar S>
    Mvec<U> operator-(const S &value, const Mvec<U> &mv){
        return (-mv) + value;
    }


    template<typename T>
    template<Scalar S>
    Mvec<T> Mvec<T>::operator-(const S &value) const {
        Mvec<T> mv(*this);
        if(value != T(0)) {
            auto it = mv.createVectorXdIfDoesNotExist(0);
            it->vec.coeffRef(0) -= value;
            mv.eraseIfNull(it);
        }
        return mv;
    }


    template<typename T>
    Mvec<T> Mvec<T>::operator^(const Mvec<T> &mv2) const {
        // Loop over non-empty grade of mv1 and mv2
        // call the right explicit unrolled outer function using the functions pointer called outerFunctionsContainer
        Mvec<T> mv3;
        for(const auto & itMv1 : this->mvData)
            for(const auto & itMv2 : mv2.mvData){
                if(itMv1.grade + itMv2.grade <= (int) algebraDimension ){
                    auto itMv3 = mv3.createVectorXdIfDoesNotExist(itMv1.grade + itMv2.grade);
                    outerFunctionsContainer<T>[itMv1.grade][itMv2.grade](itMv1.vec, itMv2.vec, itMv3->vec);
                    mv3.eraseIfNull(itMv3);
                }
            }
        return mv3;
    }

    template<typename U, Scalar S>
    Mvec<U> operator^(const S &value, const Mvec<U> &mv) {
        return  mv^value;
    }


    template<typename T>
    template<Scalar S>
    Mvec<T> Mvec<T>::operator^(const S &value) const {
        Mvec<T> mv2(*this);
        for(auto & itMv : mv2.mvData)
            itMv.vec *= T(value);
        return mv2;
    }


    template<typename T>
    Mvec<T> &Mvec<T>::operator^=(const Mvec &mv) {
        *this = *this ^ mv;
        return *this;
    }


    template<typename T>
    Mvec<T> Mvec<T>::operator|(const Mvec<T> &mv2) const{
        // Loop over non-empty grade of mv1 and mv2
        // call the right explicit unrolled inner function using the functions pointer called innerFunctionsContainer
        Mvec<T> mv3;
        for(const auto & itMv1 : this->mvData)
            for(const auto & itMv2 : mv2.mvData){

                // inner between a mv and a scalar gives 0
                if(itMv1.grade*itMv2.grade == 0)
                    continue;

                // perform the inner product
                int absGradeMv3 = std::abs((int)(itMv1.grade - itMv2.grade));
                auto itMv3 = mv3.createVectorXdIfDoesNotExist(absGradeMv3);
                innerFunctionsContainer<T>[itMv1.grade][itMv2.grade](itMv1.vec, itMv2.vec, itMv3->vec);

                // check if the result is non-zero
                if(!((itMv3->vec.array() != 0.0).any())){
                    mv3.mvData.erase(itMv3);
                    mv3.gradeBitmap &= ~(1<<absGradeMv3);
                }
            }
        return mv3;
    }


    template<typename U, Scalar S>
    Mvec<U> operator|(const S &value, const Mvec<U> &mv) {
        return Mvec<U>();
    }


    template<typename T>
    template<Scalar S>
    Mvec<T> Mvec<T>::operator|(const S &value) const {
        return Mvec<T>() ; // inner between a mv and a scalar gives 0 (empty multivector)
    }


    template<typename T>
    Mvec<T> &Mvec<T>::operator|=(const Mvec &mv) {
        *this = *this | mv;
        return *this;
    }


    template<typename T>
    Mvec<T> Mvec<T>::operator>(const Mvec<T> &mv2) const{
        // Loop over non-empty grade of mv1 and mv2
        // call the right explicit unrolled inner function using the functions pointer called innerFunctionsContainer
        Mvec<T> mv3;
        for(const auto & itMv1 : this->mvData)
            for(const auto & itMv2 : mv2.mvData){

                // right contraction constraint: gradeMv1 >= gradeMv2
                if(itMv1.grade < itMv2.grade)
                    continue;

                // perform the inner product
                int absGradeMv3 = std::abs((int)(itMv1.grade - itMv2.grade));
                auto itMv3 = mv3.createVectorXdIfDoesNotExist(absGradeMv3);
                innerFunctionsContainer<T>[itMv1.grade][itMv2.grade](itMv1.vec, itMv2.vec, itMv3->vec);

                // check if the result is non-zero
                if(!((itMv3->vec.array() != 0.0).any())){
                    mv3.mvData.erase(itMv3);
                    mv3.gradeBitmap &= ~(1<<absGradeMv3);
                }
            }

        return mv3;
    }


    template<typename U, Scalar S>
    Mvec<U> operator>(const S &value, const Mvec<U> &mv) {
        if( (mv.gradeBitmap & 1) == 0) return Mvec<U>();
        else return Mvec<U>( U(mv.mvData.begin()->vec.coeff(0) * value ) );
    }


    template<typename T>
    template<Scalar S>
    Mvec<T> Mvec<T>::operator>(const S &value) const {
        // return mv x value
        Mvec<T> mv(*this);
        for(auto & itMv : mv.mvData)
            itMv.vec *= T(value);
        return mv;
    }


    template<typename T>
    Mvec<T> Mvec<T>::operator<(const Mvec<T> &mv2) const{

        // Loop over non-empty grade of mv1 and mv2
        // call the right explicit unrolled inner function using the functions pointer called innerFunctionsContainer
        Mvec<T> mv3;
        for(const auto & itMv1 : this->mvData)
            for(const auto & itMv2 : mv2.mvData){

                // left contraction constraint: gradeMv1 <= gradeMv2
                if(itMv1.grade > itMv2.grade)
                    continue;

                // perform the inner product
                int absGradeMv3 = std::abs((int)(itMv1.grade - itMv2.grade));
                auto itMv3 = mv3.createVectorXdIfDoesNotExist(absGradeMv3);
                innerFunctionsContainer<T>[itMv1.grade][itMv2.grade](itMv1.vec, itMv2.vec, itMv3->vec);

                // check if the result is non-zero
                if(!((itMv3->vec.array() != 0.0).any())){
                    mv3.mvData.erase(itMv3);
                    mv3.gradeBitmap &= ~(1<<absGradeMv3);
                }
            }
        return mv3;
    }


    template<typename U, Scalar S>
    Mvec<U> operator<(const S &value, const Mvec<U> &mv) {
        // return mv x value
        Mvec<U> mv3(mv);
        for(auto & itMv : mv3.mvData)
            itMv.vec *= U(value);
        return mv3;
    }


    template<typename T>
    template<Scalar S>
    Mvec<T> Mvec<T>::operator<(const S &value) const {
        if( (gradeBitmap & 1) == 0) return Mvec<T>();
        else return Mvec<T>( T(mvData.begin()->vec.coeff(0) * value ) );
    }


    template<typename T>
    Mvec<T> Mvec<T>::scalarProduct(const Mvec<T> &mv2) const{
        // Loop over non-empty grade of mv1 and mv2
        // call the right explicit unrolled inner function using the functions pointer called innerFunctionsContainer
        Mvec<T> mv3;
        for(const auto & itMv1 : this->mvData)
            for(const auto & itMv2 : mv2.mvData){

                if(itMv1.grade != itMv2.grade)
                    continue;

                // perform the inner product
                int absGradeMv3 = 0;
                auto itMv3 = mv3.createVectorXdIfDoesNotExist(absGradeMv3);
                innerFunctionsContainer<T>[itMv1.grade][itMv2.grade](itMv1.vec, itMv2.vec, itMv3->vec);

                // check if the result is non-zero
                if(!((itMv3->vec.array() != 0.0).any())){
                    mv3.mvData.erase(itMv3);
                    mv3.gradeBitmap &= ~(1<<absGradeMv3);
                }
            }
        return mv3;
    }


project_singular_metric_comment_begin
    // the dual of a multivector is a cheap operation (permutation and scaling of its components), and the outer product
    // relies on the precomputed products: these "outer product with dual" operations are directly computed from their definition.
    template<typename T>
    Mvec<T> Mvec<T>::outerPrimalDual(const Mvec<T> &mv2) const{
        return (*this) ^ mv2.dual();
    }

    template<typename T>
    Mvec<T> Mvec<T>::outerDualPrimal(const Mvec<T> &mv2) const{
        return this->dual() ^ mv2;
    }

    template<typename T>
    Mvec<T> Mvec<T>::outerDualDual(const Mvec<T> &mv2) const{
        return this->dual() ^ mv2.dual();
    }

project_singular_metric_comment_end

    template<typename T>
    Mvec<T> Mvec<T>::dotProduct(const Mvec<T> &mv2) const{
        // Loop over non-empty grade of mv1 and mv2
        // call the right explicit unrolled inner function using the functions pointer called innerFunctionsContainer
        Mvec<T> mv3;
        for(const auto & itMv1 : this->mvData)
            for(const auto & itMv2 : mv2.mvData){

                // perform the inner product
                int absGradeMv3 = std::abs((int)(itMv1.grade - itMv2.grade));
                auto itMv3 = mv3.createVectorXdIfDoesNotExist(absGradeMv3);
                innerFunctionsContainer<T>[itMv1.grade][itMv2.grade](itMv1.vec, itMv2.vec, itMv3->vec);

                // check if the result is non-zero
                if(!((itMv3->vec.array() != 0.0).any())){
                    mv3.mvData.erase(itMv3);
                    mv3.gradeBitmap &= ~(1<<absGradeMv3);
                }
            }
        return mv3;
    }


    template<typename T>
    Mvec<T> Mvec<T>::operator*(const Mvec<T> &mv2) const {
        // Loop over non-empty grade of mv1 and mv2
        // call the right explicit unrolled product function using the functions pointer
        Mvec<T> mv3;
        for(const auto & itMv1 : this->mvData)
            for(const auto & itMv2 : mv2.mvData){
                project_select_recursive_geometric_product_template
                // outer product block
                unsigned int gradeOuter = itMv1.grade + itMv2.grade;
                if(gradeOuter <=  algebraDimension ){
                    auto itMv3 = mv3.createVectorXdIfDoesNotExist(gradeOuter);
                    outerFunctionsContainer<T>[itMv1.grade][itMv2.grade](itMv1.vec, itMv2.vec, itMv3->vec);
                    if(!((itMv3->vec.array() != 0.0).any())){
                        mv3.mvData.erase(itMv3);
                        mv3.gradeBitmap &= ~(1<<gradeOuter);
                    }
                }

                // inner product block
                unsigned int gradeInner = (unsigned int)std::abs((int)(itMv1.grade-itMv2.grade));
                // when the grade of one of the kvectors is zero, the inner product is the same as the outer product
                if(gradeInner != gradeOuter) {
                    auto itMv3 = mv3.createVectorXdIfDoesNotExist(gradeInner);
                    innerFunctionsContainer<T>[itMv1.grade][itMv2.grade](itMv1.vec, itMv2.vec, itMv3->vec);
                    // check if the result is non-zero
                    if(!((itMv3->vec.array() != 0.0).any())){
                        mv3.mvData.erase(itMv3);
                        mv3.gradeBitmap &= ~(1<<gradeInner);
                    }

                    // geometric product part
                    int gradeMax = (std::min)(((2*algebraDimension)-gradeOuter)+1,gradeOuter);
                    for (int gradeResult = gradeInner+2; gradeResult < gradeMax; gradeResult+=2) {
                        auto itMv3 = mv3.createVectorXdIfDoesNotExist(gradeResult);
                        geometricFunctionsContainer<T>[gradeResult][itMv1.grade][itMv2.grade](itMv1.vec, itMv2.vec, itMv3->vec);
                        // check if the result is non-zero
                        if(!((itMv3->vec.array() != 0.0).any())){
                            mv3.mvData.erase(itMv3);
                            mv3.gradeBitmap &= ~(1<<gradeResult);
                        }
                    }
                }
            }
        return mv3;
    }


    template<typename U, Scalar S>
    Mvec<U> operator*(const S &value, const Mvec<U> &mv){
        return mv * value;
    }


    template<typename T>
    template<Scalar S>
    Mvec<T> Mvec<T>::operator*(const S &value) const{
        Mvec<T> mv(*this);
        for(auto & itMv : mv.mvData)
            itMv.vec *= T(value);
        return mv;
    }


    template<typename T>
    Mvec<T> &Mvec<T>::operator*=(const Mvec &mv) {
        *this = *this * mv;
        return *this;
    }


    template<typename T>
    Mvec<T> Mvec<T>::operator/(const Mvec<T> &mv2) const {
        return *this * mv2.inv();
    }


    template<typename U, Scalar S>
    Mvec<U> operator/(const S &value, const Mvec<U> &mv){
        return value * mv.inv();
    }


    template<typename T>
    template<Scalar S>
    Mvec<T> Mvec<T>::operator/(const S &value) const {
        Mvec<T> mv(*this);
        for(auto & itMv : mv.mvData)
            itMv.vec /= value;
        return mv;
    }


    template<typename T>
    Mvec<T> &Mvec<T>::operator/=(const Mvec &mv) {
        *this = *this / mv;
        return *this;
    }

    template<typename T>
    template<Scalar S>
    Mvec<T> &Mvec<T>::operator/=(const S &value) {
        *this = *this / value;
        return *this;
    }

    template<typename T>
    Mvec<T> operator~(const Mvec<T> &mv){
        return mv.reverse();
    }


    template<typename T>
    Mvec<T> Mvec<T>::inv() const {
        T n = this->quadraticNorm();
        if(n<std::numeric_limits<T>::epsilon() && n>-std::numeric_limits<T>::epsilon())
            return Mvec<T>(); // return 0, this is was gaviewer does.
        return this->reverse() / n;
    }


    template<typename U>
    Mvec<U> operator!(const Mvec<U> &mv) {
        return mv.dual();
    }

    /// \brief defines the left contraction between two multivectors
    /// \param mv1 - a multivector
    /// \param mv2 - second operand, corresponds to a multivector
    /// \return mv1 left contraction mv2
    template<typename T>
    Mvec<T> leftContraction(const Mvec<T> &mv1, const Mvec<T> &mv2){
        return mv1 < mv2;
    }


    /// \brief defines the right contraction between two multivectors
    /// \param mv1 - a multivector
    /// \param mv2 - second operand, corresponds to a multivector
    /// \return mv1 right contraction mv2
    template<typename T>
    Mvec<T> rightContraction(const Mvec<T> &mv1, const Mvec<T> &mv2){
        return mv1 > mv2;
    }


    /// \brief returns a multivector that only contains the coefficient associated to the pseudoscalar.
    /// \return an empty Mvec if the requested element is not part of the multivector, or the multivector that contains only this element if present in the current multivector.
    template<typename T>
    Mvec<T> I() {
        Mvec<T> mvec;
        mvec[project_pseudo_scalar] = 1.0;
        return mvec;
    }

project_singular_metric_comment_begin
    /// \brief return the inverse of the pseudo scalar.
    /// \return returns a multivector corresponding to the inverse of the pseudo scalar.
    template<typename T>
    Mvec<T> Iinv() {
        Mvec<T> mvec;
        mvec[project_pseudo_scalar] = pseudoScalarInverse; // we do not return just a scalar of type T because the grade should be dimension and not 0.
        return mvec;
    }
project_singular_metric_comment_end

    // compute the dual of a multivector (i.e. mv* = mv < I, the left contraction of mv onto the pseudo-scalar)
    template<typename T>
    Mvec<T> Mvec<T>::dual() const {
        // general metric (neither diagonal nor a permutation of a diagonal matrix): the dual of a basis blade
        // is a linear combination of basis blades, it is computed from its definition
        if(!fastDualAvailable)
            return (*this) < I<T>();

        Mvec<T> mvResult;
        // for each k-vectors of the multivector
        for(auto itMv=mvData.rbegin(); itMv!=mvData.rend(); ++itMv){

            // create the dual k-vector
            Kvec<T> kvec ={itMv->vec,algebraDimension-(itMv->grade)};

            // some elements need to be permuted
            for(unsigned int i=0;i<binomialArray[itMv->grade];++i)
                kvec.vec.coeffRef(dualPermutations[itMv->grade][i]) = itMv->vec.coeff(i);

            // the inner product may involve some constant multiplucation for the dual elements
            kvec.vec = kvec.vec.cwiseProduct(dualCoefficients[itMv->grade].template cast<T>());

            // add the k-vector to the resulting multivector
            mvResult.mvData.push_back(kvec);
            mvResult.gradeBitmap |= (1 << kvec.grade);
        }
        return mvResult;
    }

    // \brief compute the reverse of a multivector
    // \return - the reverse of the multivector
    template<typename T>
    Mvec<T> Mvec<T>::reverse() const {
        Mvec<T> mv(*this);
        for(auto & itMv : mv.mvData)
            if(signReversePerGrade[itMv.grade] == -1)
                itMv.vec *= -1;
        return mv;
    }


    template<typename T>
    Mvec<T> Mvec<T>::grade(const int i) const{

        auto it = findGrade(i);
        Mvec<T> mv;

        // if not found, return the empty multivector
        if(it == mvData.end())
            return mv;

        // else return the grade 'i' data
        mv.mvData.push_back(*it);
        mv.gradeBitmap = 1 << (i);

        return mv;
    }


    template<typename T>
    std::vector<unsigned int> Mvec<T>::grades() const {
        std::vector<unsigned int> mvGrades;
        for(unsigned int i=0; i< algebraDimension+1; ++i){

            if( (gradeBitmap & (1<<i)))
                mvGrades.push_back(i);
        }

        return mvGrades;
    }


    template<typename T>
    void Mvec<T>::roundZero(const T epsilon) {
        // loop over each k-vector of the multivector
        auto itMv=mvData.begin();
        while(itMv != mvData.end()){
            // loop over each element of the k-vector
            for(unsigned int i=0; i<(unsigned int)itMv->vec.size(); ++i)
                if(fabs(itMv->vec.coeff(i)) <= epsilon)
                    itMv->vec.coeffRef(i) = 0.0;
            // if the k-vector is full of 0, remove it
            if(!((itMv->vec.array() != 0.0).any())){
                gradeBitmap = gradeBitmap - (1 << itMv->grade);
                mvData.erase(itMv++);
            }
            else ++itMv;
        }
    }


    template<typename T>
    void Mvec<T>::clear(const int grade) {

        // full erase
        if(grade < 0){
            mvData.clear();
            gradeBitmap = 0;
            return;
        }

        // partial erase
        auto iter = findGrade(grade);
        if(iter != mvData.end()) {
            gradeBitmap = gradeBitmap - (1 << iter->grade);
            iter = mvData.erase(iter);
        }
    }


    template<typename T>
    void Mvec<T>::display() const {
        if(mvData.size() == 0)
            std::cout << " null grade , null value " <<std::endl;

        for(auto itmv=mvData.begin(); itmv!=mvData.end(); ++itmv){
            std::cout << "  grade   : " << itmv->grade  << std::endl;
            std::cout << "  kvector : " << itmv->vec.transpose() << std::endl;
        }
        std::cout << std::endl;
    }


    /// \cond DEV
    template<typename U>
    void traverseKVector(std::ostream &stream, const Eigen::Matrix<U, Eigen::Dynamic, 1> kvector, unsigned int gradeMV, bool& moreThanOne ){

        // version with XOR indices
        // for the current grade, generate all the XOR indices
        std::vector<bool> booleanXorIndex(algebraDimension);
        std::fill(booleanXorIndex.begin(), booleanXorIndex.end(), false);
        // build the first xorIndex of grade 'grade' (with 'grade' values to true).
        std::fill(booleanXorIndex.begin(), booleanXorIndex.begin() + gradeMV, true);
        unsigned int positionInKVector = 0;

        do {
            // convert the vector of bool to the string containing all the basis vectors
            std::string basisBlades={};
            for(unsigned int i=0; i<algebraDimension; ++i) {
                if(booleanXorIndex[i]) {
                    basisBlades += basisVectors[i];
                }
            }

            if(kvector.coeff(positionInKVector)!= 0) {
                if(!(moreThanOne)){
                    stream<< kvector.coeff(positionInKVector) << "*e"+ basisBlades;
                    moreThanOne = true;
                }else{
                    if(kvector.coeff(positionInKVector)>0)
                        stream<< " + " << kvector.coeff(positionInKVector) << "*e" + basisBlades;
                    else
                        stream<< " - " << -kvector.coeff(positionInKVector) << "*e" + basisBlades;
                }
            }

            positionInKVector++;

        } while(std::prev_permutation(booleanXorIndex.begin(), booleanXorIndex.end())); // compute next permutation of the true values of booleanXorIndex

    }
    /// \endcond


    template<typename U>
    std::ostream &operator<<(std::ostream &stream, const Mvec<U> &mvec) {
        if(mvec.mvData.size()==0){
            stream << " 0 ";
            return stream;
        }

        bool moreThanOne = false;
        auto mvIterator = mvec.mvData.begin();

        // if the iterator corresponds to the scalar
        if(mvIterator->grade == 0){
            stream << mvIterator->vec.coeff(0);
            moreThanOne = true;
            ++mvIterator;
        }

        // for all other k-vectors of mvec
        for(; mvIterator!=mvec.mvData.end(); ++mvIterator){

            // call the function that covers a single k-vector
            traverseKVector(stream,mvIterator->vec,mvIterator->grade, moreThanOne);
        }

        if(!moreThanOne)
            stream << " 0 ";

        return stream;
    }

    // static unit multivector functions (returns a multivector with only one component)
project_static_multivector_one_component



}     /// End of Namespace

#endif // project_inclusion_guard
