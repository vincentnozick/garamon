// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// p3ga2Tools.hpp
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr

/// \file c3ga2Tools.hpp
/// \author Vincent Nozick, and Stephane Breuils
/// \brief some useful functions when using double conformal geometric algebra of R^3. Use this file if you generated the lib using the "c3ga2.conf" configuration file. Double conformal geometric algebra of R^3 is described in the following paper: Double Conformal Geometric Algebra, Robert Benjamin Easter, Eckhard Hitzer: Adv. in Appl. Clifford Algebras, 20 April 2017.



// Anti-doublon
#ifndef C3GA2_TOOLS_HPP__
#define C3GA2_TOOLS_HPP__
#pragma once

// Internal Includes
#include <c3ga2/Mvec.hpp>


/// \namespace grouping the multivectors object
namespace c3ga2{

    /// \brief build a point from a vector
    /// \param x vector component related to e1
    /// \param y vector component related to e2
    /// \param z vector component related to e3
    /// \return a multivector corresponding to a point of DPGA
    template<typename T>
    p3ga2::Mvec<T> point(const T &x, const T &y, const T &z){

        p3ga2::Mvec<T> mv;
        mv[p3ga2::E0] = 1.0;
        mv[p3ga2::E1] = x;
        mv[p3ga2::E2] = y;
        mv[p3ga2::E3] = z;

        return mv;
    }

    /// \brief build a dual point from a vector (actually, more a conjugate than a dual)
    /// \param x vector component related to e1
    /// \param y vector component related to e2
    /// \param z vector component related to e3
    /// \return a multivector corresponding to a dual point of DPGA
    template<typename T>
    p3ga2::Mvec<T> dualPoint(const T &x, const T &y, const T &z){

        p3ga2::Mvec<T> mv;
        mv[p3ga2::Ed0] = 1.0;
        mv[p3ga2::Ed1] = x;
        mv[p3ga2::Ed2] = y;
        mv[p3ga2::Ed3] = z;

        return mv;
    }

    /// \brief build a dual point from a vector (actually, more a conjugate than a dual)
    /// \param pt point of DPGA
    /// \return a multivector corresponding to a dual point of DPGA
    template<typename T>
    p3ga2::Mvec<T> dualPoint(const p3ga2::Mvec<T> &pt){

        p3ga2::Mvec<T> dualPt;
        dualPt[p3ga2::Ed0] = pt[p3ga2::E0];
        dualPt[p3ga2::Ed1] = pt[p3ga2::E1];
        dualPt[p3ga2::Ed2] = pt[p3ga2::E2];
        dualPt[p3ga2::Ed3] = pt[p3ga2::E3];

        return dualPt;
    }


    /// a quadric has the form (a b c d e f g h i j)
    /// with ax^2 + by^2 + cz^2 + dxy + exz + fyz + gxw + hyw + izw + jw^2 = 0 
    template<typename T>
    std::vector<T> ga2quadric(const p3ga2::Mvec<T> &mv){

        std::vector<T> quadric(10);
        quadric[0] = - mv[p3ga2::E0d0] / 4.0;
        quadric[1] = - mv[p3ga2::E1d1] / 4.0;
        quadric[2] = - mv[p3ga2::E2d2] / 4.0;
        quadric[9] = - mv[p3ga2::E3d3] / 4.0;

        quadric[3] = - mv[p3ga2::E0d1] / 2.0;
        quadric[4] = - mv[p3ga2::E0d2] / 2.0;
        quadric[5] = - mv[p3ga2::E2d1] / 2.0;
        quadric[6] = - mv[p3ga2::E0d3] / 2.0;
        quadric[7] = - mv[p3ga2::E3d1] / 2.0;
        quadric[8] = - mv[p3ga2::E2d3] / 2.0;

        return quadric;
    }

    /// a quadric has the form (a b c d e f g h i j)
    /// with ax^2 + by^2 + cz^2 + dxy + exz + fyz + gxw + hyw + izw + jw^2 = 0 
    template<typename T>
    p3ga2::Mvec<T> quadric2ga(const std::vector<T> &quadric){

        p3ga2::Mvec<T> mv;

        mv[p3ga2::E0d0] = - 4.0 * quadric[0];
        mv[p3ga2::E1d1] = - 4.0 * quadric[1];
        mv[p3ga2::E2d2] = - 4.0 * quadric[2];
        mv[p3ga2::E3d3] = - 4.0 * quadric[9];

        mv[p3ga2::E0d1] = mv[p3ga2::E1d0] = - 2.0 * quadric[3];
        mv[p3ga2::E0d2] = mv[p3ga2::E2d0] = - 2.0 * quadric[4];
        mv[p3ga2::E1d2] = mv[p3ga2::E2d1] = - 2.0 * quadric[5];
        mv[p3ga2::E0d3] = mv[p3ga2::E3d0] = - 2.0 * quadric[6];
        mv[p3ga2::E1d3] = mv[p3ga2::E3d1] = - 2.0 * quadric[7];
        mv[p3ga2::E2d3] = mv[p3ga2::E3d2] = - 2.0 * quadric[8];

        return mv;
    }


} // namespace

#endif // projection_inclusion_guard
