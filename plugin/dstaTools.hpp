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
        bool isCommutative= false;

        //  compare the sign of sigma_i_1 wedge sigma_k_2 with sigma_k_2 wedge sigma_i_1
        for(unsigned int i=1;i<3;++i){
            // 

        }   
        return isCommutative;
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








    /// \brief build a point from a vector
    /// \param x vector component related to e1
    /// \param y vector component related to e2
    /// \param z vector component related to e3
    /// \return a multivector corresponding to a point of DPGA
    template<typename T>
    dsta::Mvec<T> point(const T &x, const T &y, const T &z){

        dsta::Mvec<T> mv;
        mv[dsta::E0] = 1.0;
        mv[dsta::E1] = x;
        mv[dsta::E2] = y;
        mv[dsta::E3] = z;

        return mv;
    }

    /// \brief build a dual point from a vector (actually, more a conjugate than a dual)
    /// \param x vector component related to e1
    /// \param y vector component related to e2
    /// \param z vector component related to e3
    /// \return a multivector corresponding to a dual point of DPGA
    template<typename T>
    dsta::Mvec<T> dualPoint(const T &x, const T &y, const T &z){

        dsta::Mvec<T> mv;
        mv[dsta::Ed0] = 1.0;
        mv[dsta::Ed1] = x;
        mv[dsta::Ed2] = y;
        mv[dsta::Ed3] = z;

        return mv;
    }

    /// \brief build a dual point from a vector (actually, more a conjugate than a dual)
    /// \param pt point of DPGA
    /// \return a multivector corresponding to a dual point of DPGA
    template<typename T>
    dsta::Mvec<T> dualPoint(const dsta::Mvec<T> &pt){

        dsta::Mvec<T> dualPt;
        dualPt[dsta::Ed0] = pt[dsta::E0];
        dualPt[dsta::Ed1] = pt[dsta::E1];
        dualPt[dsta::Ed2] = pt[dsta::E2];
        dualPt[dsta::Ed3] = pt[dsta::E3];

        return dualPt;
    }


    /// a quadric has the form (a b c d e f g h i j)
    /// with ax^2 + by^2 + cz^2 + dxy + exz + fyz + gxw + hyw + izw + jw^2 = 0 
    template<typename T>
    std::vector<T> ga2quadric(const dsta::Mvec<T> &mv){

        std::vector<T> quadric(10);
        quadric[0] = - mv[dsta::E0d0] / 4.0;
        quadric[1] = - mv[dsta::E1d1] / 4.0;
        quadric[2] = - mv[dsta::E2d2] / 4.0;
        quadric[9] = - mv[dsta::E3d3] / 4.0;

        quadric[3] = - mv[dsta::E0d1] / 2.0;
        quadric[4] = - mv[dsta::E0d2] / 2.0;
        quadric[5] = - mv[dsta::E2d1] / 2.0;
        quadric[6] = - mv[dsta::E0d3] / 2.0;
        quadric[7] = - mv[dsta::E3d1] / 2.0;
        quadric[8] = - mv[dsta::E2d3] / 2.0;

        return quadric;
    }

    /// a quadric has the form (a b c d e f g h i j)
    /// with ax^2 + by^2 + cz^2 + dxy + exz + fyz + gxw + hyw + izw + jw^2 = 0 
    template<typename T>
    dsta::Mvec<T> quadric2ga(const std::vector<T> &quadric){

        dsta::Mvec<T> mv;

        mv[dsta::E0d0] = - 4.0 * quadric[0];
        mv[dsta::E1d1] = - 4.0 * quadric[1];
        mv[dsta::E2d2] = - 4.0 * quadric[2];
        mv[dsta::E3d3] = - 4.0 * quadric[9];

        mv[dsta::E0d1] = mv[dsta::E1d0] = - 2.0 * quadric[3];
        mv[dsta::E0d2] = mv[dsta::E2d0] = - 2.0 * quadric[4];
        mv[dsta::E1d2] = mv[dsta::E2d1] = - 2.0 * quadric[5];
        mv[dsta::E0d3] = mv[dsta::E3d0] = - 2.0 * quadric[6];
        mv[dsta::E1d3] = mv[dsta::E3d1] = - 2.0 * quadric[7];
        mv[dsta::E2d3] = mv[dsta::E3d2] = - 2.0 * quadric[8];

        return mv;
    }


} // namespace

#endif // projection_inclusion_guard
