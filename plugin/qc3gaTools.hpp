// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// qc3ga2Tools.hpp
// Authors: Vincent Nozick and Stephane Breuils 
// Contact: vincent.nozick@u-pem.fr


/// \file qc3ga2Tools.hpp
/// \author Vincent Nozick, and Stephane Breuils
/// \brief some useful functions when using Quadric Conformal Geometric Algebra of R^3. Use this file if you generated the lib using the "qc3ga.conf" configuration file. Quadric conformal geometric algebra of R^3 is described in the following paper: Stéphane Breuils, Vincent Nozick, Akihiro Sugimoto, Eckhard Hitzer, Quadric Conformal Geometric Algebra of R9,6, Advances in Applied Clifford Algebras, Springer Verlag, 2018, 28 (2), pp.35 


// Anti-doublon
#ifndef QC3GA_TOOLS_HPP__
#define QC3GA_TOOLS_HPP__
#pragma once

// External includes
#include <string>
#include <vector>

// Internal Includes
#include <qc3ga/Mvec.hpp>
#include <Eigen/Dense>        // rank in display quadric
#include <Eigen/Eigenvalues>  // eigen solver

/// \namespace grouping the multivectors object
namespace qc3ga{

    /// \brief build a point from a vector
    /// \param x vector component related to e1
    /// \param y vector component related to e2
    /// \param z vector component related to e3
    /// \return a multivector corresponding to a point of qc3ga
    template<typename T>
    qc3ga::Mvec<T> point(const T &x, const T &y, const T &z){
        qc3ga::Mvec<T> mv;
        mv[qc3ga::E1] = x;
        mv[qc3ga::E2] = y;
        mv[qc3ga::E3] = z;
        mv[qc3ga::Ei1] = 0.5*(x*x);
        mv[qc3ga::Ei2] = 0.5*(y*y);
        mv[qc3ga::Ei3] = 0.5*(z*z);
        mv[qc3ga::Ei4] = (x*y);
        mv[qc3ga::Ei5] = (x*z);
        mv[qc3ga::Ei6] = (y*z);
    		mv[qc3ga::E01] = mv[qc3ga::E02] = mv[qc3ga::E03] = 1.0;
        return mv;
    }

    template<typename T>
    qc3ga::Mvec<T> point(const Eigen::Matrix<T,3,1> &pt){
        qc3ga::Mvec<T> mv;
        mv[qc3ga::E1] = pt(0);
        mv[qc3ga::E2] = pt(1);
        mv[qc3ga::E3] = pt(2);
        mv[qc3ga::Ei1] = 0.5*(pt(0)*pt(0));
        mv[qc3ga::Ei2] = 0.5*(pt(1)*pt(1));
        mv[qc3ga::Ei3] = 0.5*(pt(2)*pt(2));
        mv[qc3ga::Ei4] = pt(0)*pt(1);
        mv[qc3ga::Ei5] = pt(0)*pt(2);
        mv[qc3ga::Ei6] = pt(1)*pt(2);
        mv[qc3ga::E01] = mv[qc3ga::E02] = mv[qc3ga::E03] = 1.0;
        return mv;
    }

    template<typename T>
    qc3ga::Mvec<T> sphere(const T &cx, const T &cy, const T &cz, const T &radius){
        qc3ga::Mvec<T> mv = qc3ga::point<T>(cx, cy, cz) - radius * radius * (qc3ga::ei1<T>() + qc3ga::ei2<T>() + qc3ga::ei3<T>()) / 6.0;
        return mv;
    }

    /// a quadric has the form (a b c d e f g h i j)
    /// with ax^2 + by^2 + cz^2 + dxy + exz + fyz + gxw + hyw + izw + jw^2 = 0 
    template<typename T>
    std::vector<T> ga2DualQuadric(const qc3ga::Mvec<T> &mv){

        std::vector<T> quadric(10);
        quadric[0] = - mv[qc3ga::E01] / 2.0;
        quadric[1] = - mv[qc3ga::E02] / 2.0;
        quadric[2] = - mv[qc3ga::E03] / 2.0;
        quadric[3] = - mv[qc3ga::E04] ;
        quadric[4] = - mv[qc3ga::E05] ;
        quadric[5] = - mv[qc3ga::E06] ;

        quadric[6] =  mv[qc3ga::E1];
        quadric[7] =  mv[qc3ga::E2];
        quadric[8] =  mv[qc3ga::E3];

        quadric[9] = - 3 * mv[qc3ga::Ei1];

        return quadric;
    }


    /// a quadric has the form (a b c d e f g h i j)
    /// with ax^2 + by^2 + cz^2 + dxy + exz + fyz + gx + hy + iz + j = 0 
    template<typename T>
    qc3ga::Mvec<T> dualQuadric2ga(const T a, const T b, const T c, const T d, const T e, const T f, const T g, const T h, const T i, const T j){
        qc3ga::Mvec<T> mv;
        mv[qc3ga::E01] = - 2.0 * a; // x^2
        mv[qc3ga::E02] = - 2.0 * b; // y^2
        mv[qc3ga::E03] = - 2.0 * c; // z^2
        mv[qc3ga::E04] = - d;       // xy
        mv[qc3ga::E05] = - e;       // xz
        mv[qc3ga::E06] = - f;       // yz
        mv[qc3ga::E1] = g;          // x
        mv[qc3ga::E2] = h;          // y
        mv[qc3ga::E3] = i;          // z

        mv[qc3ga::Ei1] = mv[qc3ga::Ei2] = mv[qc3ga::Ei3] = - j / 3.0;
        
        return mv;
    }

    /// a quadric has the form (a b c d e f g h i j)
    /// with ax^2 + by^2 + cz^2 + dxy + exz + fyz + gx + hy + iz + j = 0 
    template<typename T>
    qc3ga::Mvec<T> planeToDualGA(const T a, const T b, const T c, const T d){
        qc3ga::Mvec<T> mv;
        mv[qc3ga::E1] = a;          // x
        mv[qc3ga::E2] = b;          // y
        mv[qc3ga::E3] = c;          // z

        mv[qc3ga::Ei1] = mv[qc3ga::Ei2] = mv[qc3ga::Ei3] = - d / 3.0; // 1
        
        return mv;
    }



    /// a quadric has the form (a b c d e f g h i j)
    /// with ax^2 + by^2 + cz^2 + dxy + exz + fyz + gx + hy + iz + j = 0 
    template<typename T>
    qc3ga::Mvec<T> dualQuadric2ga(const std::vector<T> &quadric){
        return dualQuadric2ga<T>(quadric[0], quadric[1], quadric[2], quadric[3], quadric[4], quadric[5], quadric[6], quadric[7], quadric[8], quadric[9]);
    }




	
	/// a quadric has the form (a b c d e f g h i j)
    /// with ax^2 + by^2 + cz^2 + dxy + exz + fyz + gxw + hyw + izw + jw^2 = 0 
    template<typename T>
    std::vector<T> ga2Quadric(const qc3ga::Mvec<T> &mv){

        std::vector<T> quadric(10);
        quadric[0] = - mv[qc3ga::E1230102i203i304i405i506i6] / 2.0; // dual(e01)
        quadric[1] = - mv[qc3ga::E12301i10203i304i405i506i6] / 2.0; // dual(e02)
        quadric[2] = - mv[qc3ga::E12301i102i20304i405i506i6] / 2.0; // dual(e03)
        quadric[3] = - mv[qc3ga::E12301i102i203i30405i506i6]; // dual(e04)
        quadric[4] = - mv[qc3ga::E12301i102i203i304i40506i6]; // dual(e05)
        quadric[5] = - mv[qc3ga::E12301i102i203i304i405i506]; // dual(e06)
		
        quadric[6] = -mv[qc3ga::E2301i102i203i304i405i506i6]; // dual(e1)
        quadric[7] =  mv[qc3ga::E1301i102i203i304i405i506i6]; // dual(e2)
        quadric[8] = -mv[qc3ga::E1201i102i203i304i405i506i6]; // dual(e3)
		
        quadric[9] = 3 * mv[qc3ga::E123i102i203i304i405i506i6]; // dual(ei1)

        return quadric;
    }

    /// a quadric has the form (a b c d e f g h i j)
    /// with ax^2 + by^2 + cz^2 + dxy + exz + fyz + gxw + hyw + izw + jw^2 = 0
    /// defines the primal form of a conic
    template<typename T>
    qc3ga::Mvec<T> quadric2ga(const T a, const T b, const T c, const T d, const T e, const T f, const T g, const T h, const T i, const T j){
        qc3ga::Mvec<T> mv;

        mv[qc3ga::E1230102i203i304i405i506i6] = - 2.0 * a; // dual(e01) component
        mv[qc3ga::E12301i10203i304i405i506i6] = - 2.0 * b; // dual(e02) component
        mv[qc3ga::E12301i102i20304i405i506i6] = - 2.0 * c; // dual(e03) component
        mv[qc3ga::E12301i102i203i30405i506i6] = - d; // dual(e04) component
        mv[qc3ga::E12301i102i203i304i40506i6] = - e; // dual(e05) component
        mv[qc3ga::E12301i102i203i304i405i506] = - f; // dual(e06) component


        mv[qc3ga::E2301i102i203i304i405i506i6] = g; // dual(e1) component 
        mv[qc3ga::E1301i102i203i304i405i506i6] = -h; // dual(e2) component
        mv[qc3ga::E1201i102i203i304i405i506i6] = i; // dual(e3) component

        mv[qc3ga::E123i102i203i304i405i506i6] = mv[qc3ga::E12301i1i203i304i405i506i6] = mv[qc3ga::E12301i102i2i304i405i506i6] = j / 3.0; // dual(eij) component 
        
        return mv;
    }

    /// a quadric has the form (a b c d e f g h i j)
    /// with ax^2 + by^2 + cz^2 + dxy + exz + fyz + gxw + hyw + izw + jw^2 = 0
	///	defines the primal form of a conic
    template<typename T>
    qc3ga::Mvec<T> quadric2ga(const std::vector<T> &quadric){
        return qc3ga::quadric2ga<T>(quadric[0], quadric[1], quadric[2], quadric[3], quadric[4], quadric[5], quadric[6], quadric[7], quadric[8], quadric[9]);
    }


    /// \brief compute the dual of a 14-vector into a 1-vector without computing the inner product
    template<typename T>
    qc3ga::Mvec<T> primalSurfaceDualization( qc3ga::Mvec<T> &primalSurface){ // const ???
        qc3ga::Mvec<T> mv;

        mv[qc3ga::E01] = primalSurface[qc3ga::E1230102i203i304i405i506i6]; // dual(e01) component
        mv[qc3ga::E02] = primalSurface[qc3ga::E12301i10203i304i405i506i6]; // dual(e02) component
        mv[qc3ga::E03] = primalSurface[qc3ga::E12301i102i20304i405i506i6]; // dual(e03) component
        mv[qc3ga::E04] = primalSurface[qc3ga::E12301i102i203i30405i506i6]; // dual(e04) component
        mv[qc3ga::E05] = primalSurface[qc3ga::E12301i102i203i304i40506i6]; // dual(e05) component
        mv[qc3ga::E06] = primalSurface[qc3ga::E12301i102i203i304i405i506]; // dual(e06) component


        mv[qc3ga::E1]  = primalSurface[qc3ga::E2301i102i203i304i405i506i6];  // dual(e1) component 
        mv[qc3ga::E2]  = -primalSurface[qc3ga::E1301i102i203i304i405i506i6]; // dual(e2) component
        mv[qc3ga::E3]  = primalSurface[qc3ga::E1201i102i203i304i405i506i6];  // dual(e3) component

        mv[qc3ga::Ei1] = -primalSurface[qc3ga::E123i102i203i304i405i506i6];
        mv[qc3ga::Ei2] = -primalSurface[qc3ga::E12301i1i203i304i405i506i6];
        mv[qc3ga::Ei3] = -primalSurface[qc3ga::E12301i102i2i304i405i506i6];

        return mv;   
    }


    /// \brief compute the dual of a 14-vector into a 1-vector without computing the inner product
    template<typename T>
    qc3ga::Mvec<T> primalSurfaceDualization(const qc3ga::Mvec<T> &primalSurface){
        qc3ga::Mvec<T> mv;

        mv[qc3ga::E01] = primalSurface[qc3ga::E1230102i203i304i405i506i6]; // dual(e01) component
        mv[qc3ga::E02] = primalSurface[qc3ga::E12301i10203i304i405i506i6]; // dual(e02) component
        mv[qc3ga::E03] = primalSurface[qc3ga::E12301i102i20304i405i506i6]; // dual(e03) component
        mv[qc3ga::E04] = primalSurface[qc3ga::E12301i102i203i30405i506i6]; // dual(e04) component
        mv[qc3ga::E05] = primalSurface[qc3ga::E12301i102i203i304i40506i6]; // dual(e05) component
        mv[qc3ga::E06] = primalSurface[qc3ga::E12301i102i203i304i405i506]; // dual(e06) component


        mv[qc3ga::E1]  = -primalSurface[qc3ga::E2301i102i203i304i405i506i6]; // dual(e1) component 
        mv[qc3ga::E2]  = primalSurface[qc3ga::E1301i102i203i304i405i506i6]; // dual(e2) component
        mv[qc3ga::E3]  = -primalSurface[qc3ga::E1201i102i203i304i405i506i6]; // dual(e3) component

        mv[qc3ga::Ei1] = -primalSurface[qc3ga::E123i102i203i304i405i506i6];
        mv[qc3ga::Ei2] = -primalSurface[qc3ga::E12301i1i203i304i405i506i6];
        mv[qc3ga::Ei3] = -primalSurface[qc3ga::E12301i102i2i304i405i506i6];

        return mv;
    }


    /// \brief compute the normal of a surface on a point.
    /// \param surface is a 14-vector (primal form).
    /// \param point is a normalized point lying on the surface where the normal is estimated.
    /// \return a normal vector (e1,e2,e3) with L2 norm = 1 
    template<typename T>
    qc3ga::Mvec<T> surfaceNormal(qc3ga::Mvec<T> &surface, qc3ga::Mvec<T> &point){
    
      // extract the normal bivector
      qc3ga::Mvec<T> normalBivector;
      normalBivector = point ^ (primalSurfaceDualization(surface)) ;
           
      // convert the normal bivector in normal vector
      qc3ga::Mvec<T> normal;
      normal[qc3ga::E1] = -normalBivector[qc3ga::E101] - normalBivector[qc3ga::E204] - normalBivector[qc3ga::E305];  
      normal[qc3ga::E2] = -normalBivector[qc3ga::E202] - normalBivector[qc3ga::E104] - normalBivector[qc3ga::E306]; 
      normal[qc3ga::E3] = -normalBivector[qc3ga::E303] - normalBivector[qc3ga::E105] - normalBivector[qc3ga::E206]; 
 
      normal[qc3ga::Ei1] = 0.5*normal[qc3ga::E1]*normal[qc3ga::E1];  
      normal[qc3ga::Ei2] = 0.5*normal[qc3ga::E2]*normal[qc3ga::E2];  
      normal[qc3ga::Ei3] = 0.5*normal[qc3ga::E3]*normal[qc3ga::E3];

      normal[qc3ga::Ei4] = normal[qc3ga::E1]*normal[qc3ga::E2];  
      normal[qc3ga::Ei5] = normal[qc3ga::E1]*normal[qc3ga::E3];  
      normal[qc3ga::Ei6] = normal[qc3ga::E2]*normal[qc3ga::E3];  

      normal[qc3ga::E01] = normal[qc3ga::E02] = normal[qc3ga::E03] = 1.0;  
        
      return -normal;
    }

    /// \brief compute the normal of a surface on a point.
    /// \param surface is a 14-vector (primal form).
    /// \param point is a normalized point lying on the surface where the normal is estimated.
    /// \return a normal vector (e1,e2,e3), becareful the resulting normal vector is not normalized!! 
    template<typename T>
    qc3ga::Mvec<T> surfaceNormalPrimalWithoutWedge(qc3ga::Mvec<T> &surface, qc3ga::Mvec<T> &point){
    
      // extract the normal vector of the surface at the point given by point 
      // qc3ga::Mvec<T> normalBivector;
      qc3ga::Mvec<T> normal;
      normal[qc3ga::E1] = -0.5*point[qc3ga::E1]*surface[qc3ga::E01] -point[qc3ga::E2]*surface[qc3ga::E04] -point[qc3ga::E3]*surface[qc3ga::E05] + surface[qc3ga::E1];
      normal[qc3ga::E2] = -0.5*point[qc3ga::E2]*surface[qc3ga::E02] -point[qc3ga::E1]*surface[qc3ga::E04] -point[qc3ga::E3]*surface[qc3ga::E06] + surface[qc3ga::E2];
      normal[qc3ga::E3] = -0.5*point[qc3ga::E3]*surface[qc3ga::E03] -point[qc3ga::E1]*surface[qc3ga::E05] -point[qc3ga::E2]*surface[qc3ga::E06] + surface[qc3ga::E3];
        
      return -normal;
  }


 /// \brief compute the normal of a surface on a point.
    /// \param surface is a 2-vector (dual form).
    /// \param point is a normalized point lying on the surface where the normal is estimated.
    /// \return a normal vector (e1,e2,e3), becareful the resulting normal vector is not normalized!! 
    template<typename T>
    qc3ga::Mvec<T> dualSurfaceNormalWithoutWedge(qc3ga::Mvec<T> &surface, qc3ga::Mvec<T> &point){
    
      // extract the normal vector of the surface at the point given by point 
      // qc3ga::Mvec<T> normalBivector;
      qc3ga::Mvec<T> normal;
      normal[qc3ga::E1] = -point[qc3ga::E1]*surface[qc3ga::E01] -point[qc3ga::E2]*surface[qc3ga::E04] -point[qc3ga::E3]*surface[qc3ga::E05] + surface[qc3ga::E1];
      normal[qc3ga::E2] = -point[qc3ga::E2]*surface[qc3ga::E02] -point[qc3ga::E1]*surface[qc3ga::E04] -point[qc3ga::E3]*surface[qc3ga::E06] + surface[qc3ga::E2];
      normal[qc3ga::E3] = -point[qc3ga::E3]*surface[qc3ga::E03] -point[qc3ga::E1]*surface[qc3ga::E05] -point[qc3ga::E2]*surface[qc3ga::E06] + surface[qc3ga::E3];
        
      return -normal;
    }





     template<typename T>
     int quadricLineIntersection(qc3ga::Mvec<T> ptOnLine1, qc3ga::Mvec<T> ptOnLine2, qc3ga::Mvec<T> q, qc3ga::Mvec<T>& pt1, qc3ga::Mvec<T>& pt2, const T epsilon){

         // extract from the points
         qc3ga::Mvec<T> m = ptOnLine1-ptOnLine2;
         m /= std::sqrt(m[qc3ga::E1]*m[qc3ga::E1]+m[qc3ga::E2]*m[qc3ga::E2] + m[qc3ga::E3]*m[qc3ga::E3]); // Normalized support vector of the line
         ptOnLine2 = ptOnLine1 - m;

         // constants
         const T m1   = m[qc3ga::E1];
         const T m2   = m[qc3ga::E2];
         const T m3   = m[qc3ga::E3];
         const T ptl1 = ptOnLine2[qc3ga::E1];
         const T ptl2 = ptOnLine2[qc3ga::E2];
         const T ptl3 = ptOnLine2[qc3ga::E3];
         const T q1   = q[qc3ga::E1];
         const T q2   = q[qc3ga::E2];
         const T q3   = q[qc3ga::E3];
         const T q01  = q[qc3ga::E01];
         const T q02  = q[qc3ga::E02];
         const T q03  = q[qc3ga::E03];
         const T q04  = q[qc3ga::E04];
         const T q05  = q[qc3ga::E05];
         const T q06  = q[qc3ga::E06];
         const T qi1  = q[qc3ga::Ei1];
         const T qi2  = q[qc3ga::Ei2];
         const T qi3  = q[qc3ga::Ei3];



        T a = (( (-0.5*m1*m1*q01)+(-0.5*m2*m2*q02)+(-0.5*m3*m3*q03)
                + (-m1*m2*q04) + (-m1*m3*q05) +(-m2*m3*q06)));

        T c = ((-0.5*ptl1*ptl1*q01)+(-0.5*ptl2*ptl2*q02)+(-0.5*ptl3*ptl3*q03)
                +(-ptl1*ptl2*q04)+ (-ptl1*ptl3*q05) +(-ptl2*ptl3*q06) + (ptl1*q1)+ (ptl2*q2)
                + (ptl3*q3) + (-(qi1+qi2+qi3)));

        T b = ((-m1*ptl1*q01)+(-m2*ptl2*q02)+(-m3*ptl3*q03)
                + (-(m1*ptl2+m2*ptl1)*q04) + (-(m1*ptl3+m3*ptl1)*q05) 
                + (-(m2*ptl3+m3*ptl2)*q06) + (m1*q1 + m2*q2 + m3*q3));


        // With a epsilon
        if(std::fabs(a) < epsilon){
            // Flat point
            pt1 = m*(-c/b) + ptOnLine2;
            //pt2 = pt1;
            return 1;
        }

        T squaredDist = (b*b) -(4.0*a*c);

        if(squaredDist< 0.0)
            return 0; // complex solutions

        // Should do it with an epsilon
        if(std::fabs(squaredDist) < epsilon){
            // Flat point
            pt1 =(m*(-b/(2.0*(a))))+ptOnLine2;
            //pt2 = pt1;
            return 1;
        }

        pt1 =(m*((-b+std::sqrt(squaredDist))/(2.0*a)))+ptOnLine2;
        pt2 =(m*((-b-std::sqrt(squaredDist))/(2.0*a)))+ptOnLine2;

        return 2;
     }



    template<typename T>
    void quadricValueDisplay(qc3ga::Mvec<T> mv, const double epsilon){

        // remove zeros
        mv.roundZero(epsilon);

        // check unhomogeneous multivectors
        if(mv.grades().size() != 1)
            std::cout << "warning: quadricValueDisplay: not homogeneous multivector" << std::endl;

        // convert to dual quadric
        qc3ga::Mvec<T> dualQaudric;
        switch(mv.grade()) {
            case 1: // dual quadric
                dualQaudric = mv;
                break;
            case 14: // primal quadric
                dualQaudric = mv.dual();
                break;
            default: // not a quadric
                std::cout << "this multivector is not a quadric:\n" << std::endl;
                break;
        }

        // data to display
        std::string quadricComponents;    // like "ax^2 + ..."
        std::string quadricCoefficients;  // like "a=42 ..."


        // extract the quadric components
        T component;

        // a : x^2
        component = dualQaudric[qc3ga::E01] / -2.0;
        if(fabs(component) > epsilon){
            // handle the sign
            std::string sign;
            if(component < 0) sign = "- ";
            else if(quadricComponents.size() > 0) sign = "+ ";
            quadricComponents   += " + a.x^2";
            quadricCoefficients += "a = " + sign + std::to_string(fabs(component)) + "\n";
        }

        // b : y^2
        component = dualQaudric[qc3ga::E02] / -2.0;
        if(fabs(component) > epsilon){
            // handle the sign
            std::string sign;
            if(component < 0) sign = "- ";
            else if(quadricComponents.size() > 0) sign = "+ ";
            quadricComponents   += " + b.y^2";
            quadricCoefficients += "b = " + sign + std::to_string(fabs(component)) + "\n";
        }

        // c : z^2
        component = dualQaudric[qc3ga::E03] / -2.0;
        if(fabs(component) > epsilon){
            // handle the sign
            std::string sign;
            if(component < 0) sign = "- ";
            else if(quadricComponents.size() > 0) sign = "+ ";
            quadricComponents   += " + c.z^2";
            quadricCoefficients += "c = " + sign + std::to_string(fabs(component)) + "\n";
        }

        // d : xy
        component = - dualQaudric[qc3ga::E04];
        if(fabs(component) > epsilon){
            // handle the sign
            std::string sign;
            if(component < 0) sign = "- ";
            else if(quadricComponents.size() > 0) sign = "+ ";
            quadricComponents   += " + d.xy";
            quadricCoefficients += "d = " + sign + std::to_string(fabs(component)) + "\n";
        }

        // e : xz
        component = - dualQaudric[qc3ga::E05];
        if(fabs(component) > epsilon){
            // handle the sign
            std::string sign;
            if(component < 0) sign = "- ";
            else if(quadricComponents.size() > 0) sign = "+ ";
            quadricComponents   += " + e.xz";
            quadricCoefficients += "e = " + sign + std::to_string(fabs(component)) + "\n";
        }

        // f : yz
        component = - dualQaudric[qc3ga::E06];
        if(fabs(component) > epsilon){
            // handle the sign
            std::string sign;
            if(component < 0) sign = "- ";
            else if(quadricComponents.size() > 0) sign = "+ ";
            quadricComponents   += " + f.yz";
            quadricCoefficients += "f = " + sign + std::to_string(fabs(component)) + "\n";
        }

        // g : x
        component = dualQaudric[qc3ga::E1];
        if(fabs(component) > epsilon){
            // handle the sign
            std::string sign;
            if(component < 0) sign = "- ";
            else if(quadricComponents.size() > 0) sign = "+ ";
            quadricComponents   += " + g.x";
            quadricCoefficients += "g = " + sign + std::to_string(fabs(component)) + "\n";
        }

        // h : y
        component = dualQaudric[qc3ga::E2];
        if(fabs(component) > epsilon){
            // handle the sign
            std::string sign;
            if(component < 0) sign = "- ";
            else if(quadricComponents.size() > 0) sign = "+ ";
            quadricComponents   += " + h.y";
            quadricCoefficients += "h = " + sign + std::to_string(fabs(component)) + "\n";
        }

        // i : z
        component = dualQaudric[qc3ga::E3];
        if(fabs(component) > epsilon){
            // handle the sign
            std::string sign;
            if(component < 0) sign = "- ";
            else if(quadricComponents.size() > 0) sign = "+ ";
            quadricComponents   += " + i.z";
            quadricCoefficients += "i = " + sign + std::to_string(fabs(component)) + "\n";
        }

        // j : constant
        component = - (dualQaudric[qc3ga::Ei1] + dualQaudric[qc3ga::Ei2] + dualQaudric[qc3ga::Ei3]);
        if(fabs(component) > epsilon){
            // handle the sign
            std::string sign;
            if(component < 0) sign = "- ";
            else if(quadricComponents.size() > 0) sign = "+ ";
            quadricComponents   += " + j";
            quadricCoefficients += "j = " + sign + std::to_string(fabs(component)) + "\n";
        }

        // matrices delta : http://mathworld.wolfram.com/QuadraticSurface.html
        // delta 1 (rank)
        Eigen::Matrix3d delta1;
        delta1 << dualQaudric[qc3ga::E01] / -2.0, dualQaudric[qc3ga::E04] / -2.0, dualQaudric[qc3ga::E05] / -2.0,
                  dualQaudric[qc3ga::E04] / -2.0, dualQaudric[qc3ga::E02] / -2.0, dualQaudric[qc3ga::E06] / -2.0,
                  dualQaudric[qc3ga::E05] / -2.0, dualQaudric[qc3ga::E06] / -2.0, dualQaudric[qc3ga::E03] / -2.0;
        Eigen::FullPivLU<Eigen::Matrix3d> lu_decomp3(delta1);
        quadricCoefficients += "rho3    = " + std::to_string(lu_decomp3.rank()) + "\n";

        // delta 2 (rank and sign of det)
        Eigen::Matrix4d delta2;
        delta2 << dualQaudric[qc3ga::E01] / -2.0, dualQaudric[qc3ga::E04] / -2.0, dualQaudric[qc3ga::E05] / -2.0, dualQaudric[qc3ga::E1] / 2.0,
                  dualQaudric[qc3ga::E04] / -2.0, dualQaudric[qc3ga::E02] / -2.0, dualQaudric[qc3ga::E06] / -2.0, dualQaudric[qc3ga::E2] / 2.0,
                  dualQaudric[qc3ga::E05] / -2.0, dualQaudric[qc3ga::E06] / -2.0, dualQaudric[qc3ga::E03] / -2.0, dualQaudric[qc3ga::E3] / 2.0,
                  dualQaudric[qc3ga::E1]  /  2.0, dualQaudric[qc3ga::E2]  /  2.0, dualQaudric[qc3ga::E3]  /  2.0, -(dualQaudric[qc3ga::Ei1] + dualQaudric[qc3ga::Ei2] + dualQaudric[qc3ga::Ei3]) * 3.0;
        Eigen::FullPivLU<Eigen::Matrix4d> lu_decomp4(delta2);
        quadricCoefficients += "rho4    = " + std::to_string(lu_decomp4.rank()) + "\n";
        quadricCoefficients += "sign(d) = ";
        if(delta2.determinant()>0) quadricCoefficients += "+";
        if(delta2.determinant()<0) quadricCoefficients += "-";
        quadricCoefficients += "\n";

        // delta1 eigen value sign
        Eigen::EigenSolver<Eigen::Matrix3d> es(delta1,false);
        std::vector<bool> eigenValuesSignPositive;
        if(es.eigenvalues().real()(0) != 0){
            if(es.eigenvalues().real()(0) > 0) eigenValuesSignPositive.push_back(true);
            else eigenValuesSignPositive.push_back(false);
        }
        if(es.eigenvalues().real()(1) != 0){
            if(es.eigenvalues().real()(1) > 0) eigenValuesSignPositive.push_back(true);
            else eigenValuesSignPositive.push_back(false);
        }
        if(es.eigenvalues().real()(2) != 0){
            if(es.eigenvalues().real()(2) > 0) eigenValuesSignPositive.push_back(true);
            else eigenValuesSignPositive.push_back(false);
        }
        int k = 1;
        if(eigenValuesSignPositive.size() > 1)
            for(unsigned int i=1; i<eigenValuesSignPositive.size(); ++i)
                if(eigenValuesSignPositive[i] != eigenValuesSignPositive[0])
                    k = 0;

        quadricCoefficients += "k       = " + std::to_string(k) + "\n";

        // print the result
        std::cout << "\ncurrent quadric : " << quadricComponents << " = 0" << std::endl;
        std::cout << quadricCoefficients << std::endl;
    }


} // namespace

#endif // projection_inclusion_guard
