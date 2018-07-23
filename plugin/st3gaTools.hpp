// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// st3gaTools.hpp
// Authors: Vincent Nozick and Stephane Breuils
// Contact: vincent.nozick@u-pem.fr

/// \file st3gaTools.hpp
/// \author Vincent Nozick, Stephane Breuils
/// \brief some useful functions when using Space Time algebra. Use this file if you generated the lib using the "st3ga.conf" configuration file. Space Time geometric algebra is described in the following paper: 
/// Doran, Chris and Lasenby, Anthony and Gull, Stephen
/// Gravity, Gauge theories and geometric algebra, The Royal Society 1737, volume 356 (1998), pp487-582


// Anti-doublon
#ifndef ST3GA_TOOLS_HPP__
#define ST3GA_TOOLS_HPP__
#pragma once

// Internal Includes
#include <st3ga/Mvec.hpp>


/// \namespace grouping the multivectors object
namespace st3ga{

    /// \brief build the timelike bivectors basis used, ith vector
    /// note that i= 1,2,3 
    template<typename T> 
    st3ga::Mvec<T> sigma_i(const unsigned int i){
        st3ga::Mvec<T> mv;
        mv[1+1<<i]  = -1.0; // compute gamma 0 \wedge gamma i
        return mv;
    }


    /// \brief build the bivectors basis used, first timelike bivector
    template<typename T>
    st3ga::Mvec<T> sigma_1(){
        st3ga::Mvec<T> mv;
        mv[st3ga::E01] = -1.0;
        return mv;
    }


    /// \brief build the bivectors basis used, second timelike bivector 
    template<typename T>
    st3ga::Mvec<T> sigma_2(){
        st3ga::Mvec<T> mv;
        mv[st3ga::E02] = -1.0;
        return mv;
    }

    /// \brief build the bivectors basis used, third timelike bivector 
    template<typename T>
    st3ga::Mvec<T> sigma_3(){
        st3ga::Mvec<T> mv;
        mv[st3ga::E03] = -1.0;
        return mv;
    }


    /// \brief build the pseudo scalar of st3ga, using the notation of the paper
    template<typename T>
    st3ga::Mvec<T> i(){
        st3ga::Mvec<T> mv;
        mv[st3ga::E0123] = 1.0;
        return mv;
    }






    /// \brief build the projection operator associated to sigma_1 
    template<typename T>
    st3ga::Mvec<T> P_1(){
        st3ga::Mvec<T> mv;
        mv = 0.5*(1+(sigma_1<T>()));
        return mv;
    }
    
    /// \brief build the projection operator associated to sigma_2 
    template<typename T>
    st3ga::Mvec<T> P_2(){
        st3ga::Mvec<T> mv;
        mv = 0.5*(1+(sigma_2<T>()));
        return mv;
    }

    /// \brief build the projection operator associated to sigma_3 
    template<typename T>
    st3ga::Mvec<T> P_3(){
        st3ga::Mvec<T> mv;
        mv = 0.5*(1+(sigma_3<T>()));
        return mv;
    }



    /// \brief return whether the multivector mv is idempotent, i.e. mv*mv = mv
    template<typename T>
    bool isIdempotent(st3ga::Mvec<T> mv){
        return ((mv*mv)==(mv));
    }








} // namespace

#endif // projection_inclusion_guard
