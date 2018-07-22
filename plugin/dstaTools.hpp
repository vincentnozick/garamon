// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// dstaTools.hpp
// Authors: Vincent Nozick and Stephane Breuils
// Contact: vincent.nozick@u-pem.fr

/// \file dstaTools.hpp
/// \author Vincent Nozick, Stephane Breuils
/// \brief some useful functions when using double Space Time algebra. Use this file if you generated the lib using the "dsta.conf" configuration file. Double Space time geometric algebra is described in the following paper: 
/// Doran, Chris and Lasenby, Anthony and Gull, Stephen and Somaroo, Shyamal and Challinor Anthony
/// Spacetime Algebra and Electron Physics. Advances in Imaging and Electron Physics 95(2005)



// Anti-doublon
#ifndef DSTA_TOOLS_HPP__
#define DSTA_TOOLS_HPP__
#pragma once

// Internal Includes
#include <dsta/Mvec.hpp>


/// \namespace grouping the multivectors object
namespace dsta{

    /// \brief build the bivectors basis used, ith vector of the jth particle
    /// note that i= 1,2,3 and j=1,2 (first or second particle) 
    template<typename T> 
    dsta::Mvec<T> sigma_j_i(const unsigned int j, const unsigned int i){
        dsta::Mvec<T> mv;
        mv[(1<<(4*(j-1))) + (1<<((4*(j-1)) + i )) ]  = 1.0;
        return mv;
    }


    /// \brief return the coefficient associated to the four blade 
    /// note that i= 1,2,3 and j=1,2 (first or second particle) 
    template<typename T> 
    T coefficient(dsta::Mvec<T> mv, const unsigned int firstParticleIndex, const unsigned int secondParticleIndex){
        return mv[1+ (1<<4) + (1<<( firstParticleIndex)) + (1<<((4 + secondParticleIndex )))];
    }

    /// \brief build the bivectors basis used, first particle 
    template<typename T>
    dsta::Mvec<T> sigma_1_1(){
        dsta::Mvec<T> mv;
        mv[dsta::E1_0_1_1_] = 1.0;
        return mv;
    }


    /// \brief build the bivectors basis used, first particle 
    template<typename T>
    dsta::Mvec<T> sigma_1_2(){
        dsta::Mvec<T> mv;
        mv[dsta::E1_0_1_2_] = 1.0;
        return mv;
    }

    /// \brief build the bivectors basis used, first particle 
    template<typename T>
    dsta::Mvec<T> sigma_1_3(){
        dsta::Mvec<T> mv;
        mv[dsta::E1_0_1_3_] = 1.0;
        return mv;
    }


    /// \brief build the bivectors basis used, second particle 
    template<typename T>
    dsta::Mvec<T> sigma_2_1(){
        dsta::Mvec<T> mv;
        mv[dsta::E2_0_2_1_] = 1.0;
        return mv;
    }

    /// \brief build the bivectors basis used, second particle 
    template<typename T>
    dsta::Mvec<T> sigma_2_2(){
        dsta::Mvec<T> mv;
        mv[dsta::E2_0_2_2_] = 1.0;
        return mv;
    }

    /// \brief build the bivectors basis used, second particle 
    template<typename T>
    dsta::Mvec<T> sigma_2_3(){
        dsta::Mvec<T> mv;
        mv[dsta::E2_0_2_3_] = 1.0;
        return mv;
    }


    /// \brief build the pseudo scalar of the first particle, see Equation (9.6) of Section 9.1 of the reference
    template<typename T>
    dsta::Mvec<T> i_1(){
        dsta::Mvec<T> mv;
        mv[dsta::E1_0_1_1_1_2_1_3_] = 1.0;
        return mv;
    }


    /// \brief build the pseudo scalar of the second particle, see Equation (9.6) of Section 9.1 of the reference 
    template<typename T>
    dsta::Mvec<T> i_2(){
        dsta::Mvec<T> mv;
        mv[dsta::E2_0_2_1_2_2_2_3_] = 1.0;
        return mv;
    }


    // Check commutativity of the new basis, see Equation (9.5) of Section 9.1 of the reference
    template<typename T>
    bool isSigma_i_jCommutative(){

        //  compare the sign of sigma_i_1 wedge sigma_k_2 with sigma_k_2 wedge sigma_i_1
        for(unsigned int i=1;i<3;++i){ // first particle index
            for(unsigned int k=1;k<3;++i){ // second particle index
                if(i!=k){
                    dsta::Mvec<T> sigma1 = sigma_j_i<T>(1,i); // first particle blade
                    dsta::Mvec<T> sigma2 = sigma_j_i<T>(2,k); // second particle blade

                    dsta::Mvec<T> resultDirect = sigma1 ^ sigma2;
                    dsta::Mvec<T> resultPermuted = sigma2 ^ sigma1;
                    if(std::abs(coefficient(resultDirect, i, k) - coefficient(resultPermuted, i, k)) != epsilon) 
                        return false;
                }
            }
        }   
        return true;
    }


    /// \brief build the projection operator which is called E, see Equation (9.10) of Section 9.1 of the reference
    template<typename T>
    dsta::Mvec<T> E(){
        dsta::Mvec<T> mv;
        mv = 0.5*(1-(i_1<T>()*sigma_1_3<T>() + i_2<T>()*sigma_2_3<T>()));
        return mv;
    }
    

    /// \brief build the projection operator which is called E, see Equation (9.12) of Section 9.1 of the reference
    template<typename T>
    dsta::Mvec<T> J(){
        dsta::Mvec<T> mv;
        mv = (E<T>()*i_1<T>()*sigma_1_3<T>());
        return mv;
    }







} // namespace

#endif // projection_inclusion_guard
