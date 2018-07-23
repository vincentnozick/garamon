// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// st3ga2Tools.hpp
// Authors: Vincent Nozick and Stephane Breuils
// Contact: vincent.nozick@u-pem.fr

/// \file st3ga2Tools.hpp
/// \author Vincent Nozick, Stephane Breuils
/// \brief some useful functions when using double Space Time algebra. Use this file if you generated the lib using the "st3ga2.conf" configuration file. Double Space time geometric algebra is described in the following paper: 
/// Doran, Chris and Lasenby, Anthony and Gull, Stephen and Somaroo, Shyamal and Challinor Anthony
/// Spacetime Algebra and Electron Physics. Advances in Imaging and Electron Physics 95(2005)



// Anti-doublon
#ifndef ST3GA2_TOOLS_HPP__
#define ST3GA2_TOOLS_HPP__
#pragma once

// Internal Includes
#include <st3ga2/Mvec.hpp>


/// \namespace grouping the multivectors object
namespace st3ga2{

    /// \brief build the bivectors basis used, ith vector of the jth particle
    /// note that i= 1,2,3 and j=1,2 (first or second particle) 
    template<typename T> 
    st3ga2::Mvec<T> sigma_j_i(const unsigned int j, const unsigned int i){
        st3ga2::Mvec<T> mv;
        mv[(1<<(4*(j-1))) + (1<<((4*(j-1)) + i )) ]  = -1.0;
        return mv;
    }


    /// \brief return the coefficient associated to the four blade 
    /// note that i= 1,2,3 and j=1,2 (first or second particle) 
    template<typename T> 
    T coefficient(st3ga2::Mvec<T> mv, const unsigned int firstParticleIndex, const unsigned int secondParticleIndex){
        return mv[1+ (1<<4) + (1<<( firstParticleIndex)) + (1<<((4 + secondParticleIndex )))];
    }

    /// \brief build the bivectors basis used, first particle 
    template<typename T>
    st3ga2::Mvec<T> sigma_1_1(){
        st3ga2::Mvec<T> mv;
        mv[st3ga2::E10_11_] = -1.0;
        return mv;
    }


    /// \brief build the bivectors basis used, first particle 
    template<typename T>
    st3ga2::Mvec<T> sigma_1_2(){
        st3ga2::Mvec<T> mv;
        mv[st3ga2::E10_12_] = -1.0;
        return mv;
    }

    /// \brief build the bivectors basis used, first particle 
    template<typename T>
    st3ga2::Mvec<T> sigma_1_3(){
        st3ga2::Mvec<T> mv;
        mv[st3ga2::E10_13_] = -1.0;
        return mv;
    }


    /// \brief build the bivectors basis used, second particle 
    template<typename T>
    st3ga2::Mvec<T> sigma_2_1(){
        st3ga2::Mvec<T> mv;
        mv[st3ga2::E20_21_] = -1.0;
        return mv;
    }

    /// \brief build the bivectors basis used, second particle 
    template<typename T>
    st3ga2::Mvec<T> sigma_2_2(){
        st3ga2::Mvec<T> mv;
        mv[st3ga2::E20_22_] = -1.0;
        return mv;
    }

    /// \brief build the bivectors basis used, second particle 
    template<typename T>
    st3ga2::Mvec<T> sigma_2_3(){
        st3ga2::Mvec<T> mv;
        mv[st3ga2::E20_23_] = -1.0;
        return mv;
    }


    /// \brief build the pseudo scalar of the first particle, see Equation (9.6) of Section 9.1 of the reference
    template<typename T>
    st3ga2::Mvec<T> i_1(){
        st3ga2::Mvec<T> mv;
        mv[st3ga2::E10_11_12_13_] = 1.0;
        return mv;
    }


    /// \brief build the pseudo scalar of the second particle, see Equation (9.6) of Section 9.1 of the reference 
    template<typename T>
    st3ga2::Mvec<T> i_2(){
        st3ga2::Mvec<T> mv;
        mv[st3ga2::E20_21_22_23_] = 1.0;
        return mv;
    }


    // Check commutativity of the new basis, see Equation (9.5) of Section 9.1 of the reference
    template<typename T>
    bool isSigma_i_jCommutative(const T &epsilon = 1.0e-7){

        //  compare the sign of sigma_i_1 wedge sigma_k_2 with sigma_k_2 wedge sigma_i_1
        for(unsigned int i=1;i<3;++i){ // first particle index
            for(unsigned int k=1;k<3;++k){ // second particle index
                if(i!=k){
                    st3ga2::Mvec<T> sigma1 = sigma_j_i<T>(1,i); // first particle blade
                    st3ga2::Mvec<T> sigma2 = sigma_j_i<T>(2,k); // second particle blade

                    st3ga2::Mvec<T> resultDirect = sigma1 ^ sigma2;
                    st3ga2::Mvec<T> resultPermuted = sigma2 ^ sigma1;
                    if(std::abs(coefficient(resultDirect, i, k) - coefficient(resultPermuted, i, k)) > epsilon) 
                        return false;
                }
            }
        }   
        return true;
    }




    /// \brief build the projection operator which is called E, see Equation (9.10) of Section 9.1 of the reference
    template<typename T>
    st3ga2::Mvec<T> E(){
        st3ga2::Mvec<T> mv;
        mv = 0.5*(1-(i_1<T>()*sigma_1_3<T>() * i_2<T>()*sigma_2_3<T>()));
        return mv;
    }
    
    /// \brief return whether the multivector mv is idempotent, i.e. mv*mv = mv
    template<typename T>
    bool isIdempotent(st3ga2::Mvec<T> mv){
        return ((mv*mv)==(mv));
    }


    /// \brief build the projection operator which is called E, see Equation (9.12) of Section 9.1 of the reference
    template<typename T>
    st3ga2::Mvec<T> J(){
        st3ga2::Mvec<T> mv;
        mv = (E<T>()*i_1<T>()*sigma_1_3<T>());
        return mv;
    }







} // namespace

#endif // projection_inclusion_guard
