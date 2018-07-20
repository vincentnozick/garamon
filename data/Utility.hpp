// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// Utility.hpp
// This file is part of the Garamon for project_namespace.
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

/// \file Utility.hpp
/// \author Stephane Breuils, Vincent Nozick
/// \brief basic and generic tools .


#ifndef project_inclusion_guard
#define project_inclusion_guard
#pragma once

/*!
 * @namespace project_namespace
 */
namespace project_namespace{


    template<typename T, int N>
    struct ArrayIndices{
        T ind[N];
        int sign;
        constexpr T& operator[](int i){return ind[i];}
        constexpr const T& operator[](int i) const {return ind[i];}
    };

    template<int ... paramIdx>
    constexpr auto sortAndDetermineSign(){
        // recompute the grade useful for the sequel
        constexpr int gd = sizeof...(paramIdx);

        ArrayIndices<int, gd> a = {paramIdx...};


        for (int i = 0;  i < gd - 1;  i++) {
            for (int j = 0; j < gd - i - 1; j++) {
                if (a.ind[j] > a.ind[j+1])
                {
                    a.sign *= -1;
                    int temp  = a[j];
                    a[j] = a[j+1];
                    a[j+1]= temp;
                }
            }
        }
        return a;
    }


    /*!
    * Compute the binomial coefficient
    * @tparam n
    * @tparam k
    */
    template<int n, int k>
    struct Binomial
    {
        const static int value =  (Binomial<n-1,k-1>::value + Binomial<n-1,k>::value);
    };

    template<>
    struct Binomial<0,0>
    {
        const static int value = 1;
    };

    template<int n>
    struct Binomial<n,0>
    {
        const static int value = 1;
    };

    template<int n>
    struct Binomial<0,n>
    {
        const static int value = 0;
    };

    template<int n>
    struct Binomial<n,n>
    {
        const static int value = 1;
    };

    template<const int n, const int k>
    constexpr const static int binomial()
    {
        return Binomial<n,k>::value;
    }

    constexpr unsigned int factorial(int n)
    {
        return (n <= 1) ? 1 : n * factorial(n - 1);
    }

    constexpr unsigned int bin_coeff(int n, int k)
    {
        return (k>n)?0:factorial(n) / factorial(n - k) / factorial(k);
    }


   /*!
    * From a list of vectors, compute the index in a homogeneous k-vectors
    * @tparam T - same type as Arguments, i.e. integer
    * @tparam Arguments - type Variadic
    * @param firstIdx - the first index
    * @param a - variadic list of indices
    * @example - from computeIdxVariadic(1,2,4) = 1
    * @return the index in the VectorXd corresponding to the variadic list of index.
    */
    template<int dimens, int grad, int first, int... Arguments>
    struct idxVariadic{
        const static int value = Binomial<dimens-first,grad>::value + idxVariadic<dimens,grad-1,Arguments...>::value;
    };

  /*!
   * specialized version of computeIdxVariadic;
   * @tparam T - corresponds to the type of index
   * @param last - the last index of the basis blade
   * @return the index in the VectorXd corresponding to the variadic list of index.
   */
    template<int dimens, int grad, int last>
    struct idxVariadic<dimens,grad,last>{
        const static int value = 0;
    };


    template<typename T>
    constexpr int computeIdxFromList(const int dimens, int grad){
        return 0;
    }


    template<typename T>
    constexpr int computeIdxFromList(const int dimens, int grad, T last){
        return bin_coeff(dimens-last,grad);
    }

    template<typename T, typename... ListBlades>
    constexpr int computeIdxFromList(const int dimens, int grad, T first, ListBlades ...list){
        return bin_coeff(dimens-first,grad) + computeIdxFromList(dimens,grad-1,list...);
    }


    /*
     * End of Namespace
     */
}



#endif // project_inclusion_guard
