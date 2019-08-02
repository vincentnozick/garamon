// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// cgaTools.hpp
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr

/// \file c3gaTools.hpp
/// \author Stephane Breuils, Vincent Nozick
/// \brief some useful functions when using conformal geometric algebra of R^3. Use this file if you generated the lib using the "c3ga.conf" configuration file.


// Anti-doublon
#ifndef C3GA_TOOLS_HPP__
#define C3GA_TOOLS_HPP__
#pragma once

// Internal Includes
#include <c3ga/Mvec.hpp>


/// \namespace grouping the multivectors object
namespace c3ga{

    /// \brief build a point from a vector
    /// \param x vector component related to e1
    /// \param y vector component related to e2
    /// \param z vector component related to e3
    /// \return a multivector corresponding to a point p = e0 +  x e1 + y e2 + z e3 + 0.5 || (x e1 + y e2 + z e3) ||^2 einf
    template<typename T>
    c3ga::Mvec<T> point(const T &x, const T &y, const T &z){

        c3ga::Mvec<T> mv;
        mv[c3ga::E1] = x;
        mv[c3ga::E2] = y;
        mv[c3ga::E3] = z;
        mv[c3ga::Ei] = 0.5 * mv.quadraticNorm();
        mv[c3ga::E0] = 1.0;

        return mv;
    }

    /// \brief build a point from a vector
    /// \param vec is should be a multivector of grade 1 vec = v1 e1 + v2 e2 + v3 e3. If vec has other components, they will be ignored during the execution.
    /// \return a multivector corresponding to a point p = e0 + v1 e1 + v2 e2 + v3 e3 + 0.5 || vec ||^2 einf
    template<typename T>
    c3ga::Mvec<T> point(const c3ga::Mvec<T> &vec){
        return point(vec[c3ga::E1], vec[c3ga::E2], vec[c3ga::E3]);
    }

    /// \brief build a dual sphere from a center and a radius
    /// \param centerX dual sphere center component related to e1
    /// \param centerY dual sphere center component related to e2
    /// \param centerZ dual sphere center component related to e3
    /// \param radius of the sphere
    /// \return a multivector corresponding to a dual sphere s = center - 0.5 radius ei, with center being e0 +  x e1 + y e2 + z e3 + 0.5 || (x e1 + y e2 + z e3) ||^2 einf.
    template<typename T>
    c3ga::Mvec<T> dualSphere(const T &centerX, const T &centerY, const T &centerZ, const T &radius){
        c3ga::Mvec<T> dualSphere = point(centerX,centerY,centerZ);
        dualSphere[c3ga::Ei] -= 0.5*radius;
        return dualSphere;
    }

    /// \brief extract the center and radius of a dual sphere
    /// \param dualSphere the dual sphere
    /// \param radius positive number for real spheres and negative for imaginary spheres.
    /// \param center center of the dual sphere
    template<typename T>
    void radiusAndCenterFromDualSphere(const c3ga::Mvec<T> &dualSphere, T &radius, c3ga::Mvec<T> &center){
        radius = dualSphere | dualSphere;
        center = dualSphere;
        dualSphere[c3ga::Ei] += 0.5*radius*radius;
    }

    /// \brief interpret the nature of the geometric object (line, circle, pair point, ...)
    /// \param multivector: the multivector to be studied
    /// \todo ... todo :)
    template<typename T>
    void whoAmI(const c3ga::Mvec<T> &mv){

        std::vector<unsigned int> grades_ = mv.grades();
        if(grades_.size() == 0) {
            std::cout << "null vector" << std::endl;
            return;
        }

        if(grades_.size() == 1)
            std::cout << grades_[0] << "-vector (homogeneous)" << std::endl;
        else
            std::cout << "non-homogeous multivector" << std::endl;
    }


    /// \brief extract 2 points pt1 and pt2 from a pair point p = pt1 ^ pt2
    /// \param pairPoint implicitly contains 2 points
    /// \param epsilon is the minimum threshold to specify if 2 points are disjoint
    /// \return a list of 2 points (if they are disjoint) or a single point.
    template<typename T>
    std::vector<c3ga::Mvec<T>> extractPairPoint(const c3ga::Mvec<T> &pairPoint, const T &epsilon = 1.0e-7){

        std::vector<c3ga::Mvec<T>> points;
        T innerSqrt = sqrt(pairPoint | pairPoint);
        if(innerSqrt < epsilon)
            points.push_back(pairPoint / pairPoint[c3ga::E0]);
        else {
            points.push_back( (pairPoint+innerSqrt)/ pairPoint[c3ga::E0]);
            points.push_back( (pairPoint-innerSqrt)/ pairPoint[c3ga::E0]);
        }
        return points;
    }


	/// \brief compute the normal of a surface on a point.
	/// \param surface is a quad vector (sphere of plane).
	/// \param point is a normalized point (e0 = 1) lying on the surface where the normal is estimated.
	/// \return a normal vector (e1,e2,e3) with L2 norm = 1 
    template<typename T>
    c3ga::Mvec<T> surfaceNormal(c3ga::Mvec<T> &surface, c3ga::Mvec<T> &point){
	    cga::Mvec<T> normal;
	    normal[c3ga::E1] = - point[c3ga::E1] * surface[c3ga::E0123] / point[c3ga::E0] + surface[c3ga::E023i];  
	    normal[c3ga::E2] = - point[c3ga::E2] * surface[c3ga::E0123] / point[c3ga::E0] - surface[c3ga::E013i]; 
	    normal[c3ga::E3] = - point[c3ga::E3] * surface[c3ga::E0123] / point[c3ga::E0] + surface[c3ga::E012i]; 

	    normal = normal / (double) sqrt(normal[c3ga::E1]*normal[c3ga::E1] + normal[c3ga::E2]*normal[c3ga::E2] + normal[c3ga::E3]*normal[c3ga::E3]);

	    return normal;
    }


} // namespace

#endif // projection_inclusion_guard
