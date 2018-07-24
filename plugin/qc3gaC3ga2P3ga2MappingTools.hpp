// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// p3ga2Tools.hpp
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr

/// \file qc3gaC3ga2P3ga2MappingTools.hpp
/// \author Stephane Breuils, Vincent Nozick
/// \brief some useful functions when using the mapping between double projective geometric algebra of R^3 (c3ga2).
///  Use this file if you generated the lib using the "p3ga2.conf" configuration file. Double projective geometric algebra of R^3 is described in the following paper: Juan, D., Goldman, R., Mann, S.: Modeling 3D geometry in the Clifford Algebra R 4,4 . Adv. Appl. Clifford Algebras 27(4), 3039–3062 (2017)



// Anti-doublon
#ifndef QC3GAC3GA2P3GA2MAPPINGTOOLS_HPP__
#define QC3GAC3GA2P3GA2MAPPINGTOOLS_HPP__
#pragma once

// Internal Includes
#include <p3ga2/Mvec.hpp> // double  projective geometric algebra
#include <c3ga2/Mvec.hpp> // double  conformal  geometric algebra
#include <qc3ga/Mvec.hpp> // quadric conformal  geometric algebra





    //////////////////////////////////////// P3GA2QUADRIC //////////////////////////////////////////
namespace p3ga2{
    /// \brief build a point from a vector
    /// \param x vector component related to e1
    /// \param y vector component related to e2
    /// \param z vector component related to e3
    /// \return a multivector corresponding to a point of DPGA
    template<typename T>
    p3ga2::Mvec<T> point(const T &x, const T &y, const T &z){

        p3ga2::Mvec<T> mv;
        mv[p3ga2::E0] = x;
        mv[p3ga2::E1] = y;
        mv[p3ga2::E2] = z;
        mv[p3ga2::E3] = 1.0;

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
        mv[p3ga2::Ed0] = x;
        mv[p3ga2::Ed1] = y;
        mv[p3ga2::Ed2] = z;
        mv[p3ga2::Ed3] = 1.0;

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



//////////////////////////////////////// C3GA2 QUADRIC //////////////////////////////////////////
namespace c3ga2{

   /// \brief build a point from a vector, see Equation (5.25) of the considered paper
    /// \param x vector component 
    /// \param y vector component 
    /// \param z vector component 
    /// \return a multivector corresponding to a point of DCGA
    template<typename T>
    c3ga2::Mvec<T> point(const T &x, const T &y, const T &z){
       c3ga2::Mvec<T> c3ga2Point;
        
       c3ga2::Mvec<T> mvCGA1;
        c3ga2::Mvec<T> mvCGA2;
        mvCGA1[c3ga2::E01] = mvCGA2[c3ga2::E02] = 1.0;
        mvCGA1[c3ga2::E1]  = mvCGA2[c3ga2::E4]  = x;
        mvCGA1[c3ga2::E2]  = mvCGA2[c3ga2::E5]  = y;
        mvCGA1[c3ga2::E3]  = mvCGA2[c3ga2::E6]  = z;
        mvCGA1[c3ga2::Ei1] = mvCGA2[c3ga2::Ei2] = 0.5*(x*x+y*y+z*z);

       c3ga2Point = mvCGA1^mvCGA2;
        
       return c3ga2Point;
    }


        /// \brief build Txsquared
    /// \return the multivector corresponding to the Txsquared
    template<typename T>
    c3ga2::Mvec<T> Tx2(){
       c3ga2::Mvec<T> mv;
       mv[c3ga2::E14] = -1.0;
        
       return mv;
    }

    /// \brief build Tysquared
    /// \return the multivector corresponding to the Tysquared
    template<typename T>
    c3ga2::Mvec<T> Ty2(){
        c3ga2::Mvec<T> mv;
        mv[c3ga2::E25] = -1.0;
            
        return mv;
    }

    /// \brief build Tzsquared
    /// \return the multivector corresponding to the Tzsquared
    template<typename T>
    c3ga2::Mvec<T> Tz2(){
        c3ga2::Mvec<T> mv;
        mv[c3ga2::E36] = -1.0;
        return mv;
    }


    /// \brief build the cross term operator Txy
    /// \return the multivector corresponding to the operator Txy
    template<typename T>
    c3ga2::Mvec<T> Txy(){
        c3ga2::Mvec<T> mv;
        mv[c3ga2::E15] = -0.5;
        mv[c3ga2::E24] = -0.5;
        return mv;
    }


    /// \brief build the cross term operator Tzx
    /// \return the multivector corresponding to the operator Tzx
    template<typename T>
    c3ga2::Mvec<T> Tzx(){
        c3ga2::Mvec<T> mv;
        mv[c3ga2::E16] = -0.5;
        mv[c3ga2::E34] = -0.5;
        return mv;
    }

    /// \brief build the cross term operator Tyz
    /// \return the multivector corresponding to the operator Tyz
    template<typename T>
    c3ga2::Mvec<T> Tyz(){
        c3ga2::Mvec<T> mv;
        mv[c3ga2::E35] = -0.5;
        mv[c3ga2::E26] = -0.5;
        return mv;
    }


    /// \brief build the single term operator Tx
    /// \return the multivector corresponding to the operator Tx
    template<typename T>
    c3ga2::Mvec<T> Tx(){
        c3ga2::Mvec<T> mv;
        mv[c3ga2::E1i2] = 0.5;
        mv[c3ga2::Ei14] = 0.5;
        return mv;
    }


    /// \brief build the single term operator Ty
    /// \return the multivector corresponding to the operator Ty
    template<typename T>
    c3ga2::Mvec<T> Ty(){
        c3ga2::Mvec<T> mv;
        mv[c3ga2::E2i2] = 0.5;
        mv[c3ga2::Ei15] = 0.5;
        return mv;
    }

    /// \brief build the single term operator Tz
    /// \return the multivector corresponding to the operator Tz
    template<typename T>
    c3ga2::Mvec<T> Tz(){
        c3ga2::Mvec<T> mv;
        mv[c3ga2::E3i2] = 0.5;
        mv[c3ga2::Ei16] = 0.5;
        return mv;
    }

    /// \brief build the unit term operator Tz
    /// \return the multivector corresponding to the operator T1
    template<typename T>
    c3ga2::Mvec<T> T1(){
        c3ga2::Mvec<T> mv;
        mv[c3ga2::Ei1i2] = -1.0;
        return mv;
    }



    /// \brief build the reciprocal squared term operator Tx2_reciprocal
    /// \return the multivector corresponding to the reciprocal extraction operator Tx2
    template<typename T>
    c3ga2::Mvec<T> Tx2_reciprocal(){
        return -Tx2<T>();
    }

    /// \brief build the reciprocal squared term operator Ty2_reciprocal
    /// \return the multivector corresponding to the reciprocal extraction operator Ty2
    template<typename T>
    c3ga2::Mvec<T> Ty2_reciprocal(){
        return -Ty2<T>();
    }

    /// \brief build the reciprocal squared term operator Tz2_reciprocal
    /// \return the multivector corresponding to the reciprocal extraction operator Tz2
    template<typename T>
    c3ga2::Mvec<T> Tz2_reciprocal(){
        return -Tz2<T>();
    }

    /// \brief build the reciprocal cross term operator Txy_reciprocal
    /// \return the multivector corresponding to the reciprocal extraction operator Txy
    template<typename T>
    c3ga2::Mvec<T> Txy_reciprocal(){
        return -2*Txy<T>();
    }

    /// \brief build the reciprocal cross term operator Tzx_reciprocal
    /// \return the multivector corresponding to the reciprocal extraction operator Tzx
    template<typename T>
    c3ga2::Mvec<T> Tzx_reciprocal(){
        return -2*Tzx<T>();
    }

    /// \brief build the reciprocal cross term operator Tzx_reciprocal
    /// \return the multivector corresponding to the reciprocal extraction operator Tzx
    template<typename T>
    c3ga2::Mvec<T> Tyz_reciprocal(){
        return -2*Tyz<T>();
    }

    /// \brief build the reciprocal single term operator Tx_reciprocal
    /// \return the multivector corresponding to the reciprocal extraction operator Tx
    template<typename T>
    c3ga2::Mvec<T> Tx_reciprocal(){
        c3ga2::Mvec<T> mv;
        mv[c3ga2::E102] = 1.0;
        mv[c3ga2::E014] = 1.0;
        return mv;
    }

    /// \brief build the reciprocal single term operator Ty_reciprocal
    /// \return the multivector corresponding to the reciprocal extraction operator Ty
    template<typename T>
    c3ga2::Mvec<T> Ty_reciprocal(){
        c3ga2::Mvec<T> mv;
        mv[c3ga2::E202] = 1.0;
        mv[c3ga2::E015] = 1.0;
        return mv;
    }

    /// \brief build the reciprocal single term operator Tz_reciprocal
    /// \return the multivector corresponding to the reciprocal extraction operator Tz
    template<typename T>
    c3ga2::Mvec<T> Tz_reciprocal(){
        c3ga2::Mvec<T> mv;
        mv[c3ga2::E302] = 1.0;
        mv[c3ga2::E016] = 1.0;
        return mv;
    }

    /// \brief build the reciprocal unit term operator Tz
    /// \return the multivector corresponding to the reciprocal extraction operator T1
    template<typename T>
    c3ga2::Mvec<T> T1_reciprocal(){
        c3ga2::Mvec<T> mv;
        mv[c3ga2::E0102] = 1.0;
        return mv;
    }



    /// \brief build the cyclide component extraction operator
    /// \return the multivector corresponding to the operator Tt4
    template<typename T>
    c3ga2::Mvec<T> Tt4(){
        c3ga2::Mvec<T> mv;
        mv[c3ga2::E0102] = -4.0;
        return mv;
    }

    /// \brief build the cyclide component extraction operator
    /// \return the multivector corresponding to the operator Tt2
    template<typename T>
    c3ga2::Mvec<T> Tt2(){
        c3ga2::Mvec<T> mv;
        mv[c3ga2::Ei102] = -1.0;
        mv[c3ga2::E01i2] = -1.0;
        return mv;
    }

    /// \brief build the cyclide component extraction operator
    /// \return the multivector corresponding to the operator Txt2
    template<typename T>
    c3ga2::Mvec<T> Txt2(){
        c3ga2::Mvec<T> mv;
        mv[c3ga2::E102] = 1.0;
        mv[c3ga2::E014] = 1.0;
        return mv;
    }

    /// \brief build the cyclide component extraction operator
    /// \return the multivector corresponding to the operator Tyt2
    template<typename T>
    c3ga2::Mvec<T> Tyt2(){
        c3ga2::Mvec<T> mv;
        mv[c3ga2::E202] = 1.0;
        mv[c3ga2::E015] = 1.0;
        return mv;
    }

    /// \brief build the cyclide component extraction operator
    /// \return the multivector corresponding to the operator Tzt2
    template<typename T>
    c3ga2::Mvec<T> Tzt2(){
        c3ga2::Mvec<T> mv;
        mv[c3ga2::E302] = 1.0;
        mv[c3ga2::E016] = 1.0;
        return mv;
    }





    /// a quadric has the form (a b c d e f g h i j)
    /// with ax^2 + by^2 + cz^2 + dxy + exz + fyz + gxw + hyw + izw + jw^2 = 0
    /// \brief construct the multivector representing the quadric whose coefficients are {a,b,c,d,e,f,g,h,i,j}
    template<typename T>
    c3ga2::Mvec<T> quadric2ga(const std::vector<T>& quadric){
        c3ga2::Mvec<T> mv = quadric[0]*Tx2<T>() + quadric[1]*Ty2<T>() + quadric[2]*Tz2<T>() + quadric[3]*Txy<T>() + quadric[4]*Tzx<T>() + quadric[5]*Tyz<T>()  
              + quadric[6]*Tx<T>()  + quadric[7]*Ty<T>()  + quadric[8]*Tz<T>()  + quadric[9]*T1<T>();
    
        return mv;
    }


    /// a quadric has the form (a b c d e f g h i j)
    /// with ax^2 + by^2 + cz^2 + dxy + exz + fyz + gxw + hyw + izw + jw^2 = 0
    /// \brief construct the quadric whose bivector multivector is mv, to achieve that we use the reciprocal bivector operators  
    template<typename T>
    std::vector<T> ga2quadric(c3ga2::Mvec<T> mv ){
        std::vector<T> quadric(10);
        quadric[0] = (Tx2_reciprocal<T>()|mv); // a coefficient
        quadric[1] = (Ty2_reciprocal<T>()|mv); // b coefficient
        quadric[2] = (Tz2_reciprocal<T>()|mv); // c coefficient
        quadric[3] = (Txy_reciprocal<T>()|mv); // d coefficient
        quadric[4] = (Tzx_reciprocal<T>()|mv); // e coefficient
        quadric[5] = (Tyz_reciprocal<T>()|mv); // f coefficient
        quadric[6] = (Tx_reciprocal<T>()|mv);  // g coefficient
        quadric[7] = (Ty_reciprocal<T>()|mv);  // h coefficient
        quadric[8] = (Tz_reciprocal<T>()|mv);  // i coefficient
        quadric[9] = (T1_reciprocal<T>()|mv);  // j coefficient
        return quadric;
    }

} // end namespace

namespace qc3ga{
    //////////////////////////////////////// QC3GA QUADRIC //////////////////////////////////////////
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

        mv[qc3ga::E01] = - 2.0 * a;
        mv[qc3ga::E02] = - 2.0 * b;
        mv[qc3ga::E03] = - 2.0 * c;
        mv[qc3ga::E04] = - d;
        mv[qc3ga::E05] = - e;
        mv[qc3ga::E06] = - f;

        mv[qc3ga::E1] = g;
        mv[qc3ga::E2] = h;
        mv[qc3ga::E3] = i;

        mv[qc3ga::Ei1] = mv[qc3ga::Ei2] = mv[qc3ga::Ei3] = - j / 3.0;
        
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
    std::vector<T> ga2quadric(const qc3ga::Mvec<T> &mv){

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
    /// defines the primal form of a conic
    template<typename T>
    qc3ga::Mvec<T> quadric2ga(const std::vector<T> &quadric){
        return qc3ga::quadric2ga<T>(quadric[0], quadric[1], quadric[2], quadric[3], quadric[4], quadric[5], quadric[6], quadric[7], quadric[8], quadric[9]);
    }

} // end namespace


    /////////////////////////// MAPPING //////////////////////////

    /// \brief transform a c3ga2 quadric to a qc3ga one
    template<typename T>
    qc3ga::Mvec<T> c3ga2quadric2qc3ga(const c3ga2::Mvec<T> & mvQuadric){
    std::vector<T> quadCoeff = c3ga2::ga2quadric<T>(mvQuadric);
        return qc3ga::dualQuadric2ga<T>(quadCoeff);
    }

    /// \brief transform a qc3ga quadric to a c3ga2 one
    template<typename T>
    c3ga2::Mvec<T> qc3gaquadric2c3ga2(const qc3ga::Mvec<T> & mvQuadric){
    std::vector<T> quadCoeff = qc3ga::ga2DualQuadric<T>(mvQuadric);
        return c3ga2::quadric2ga<T>(quadCoeff);
    }


    /// \brief transform a p3ga2 quadric to a qc3ga one
    template<typename T>
    qc3ga::Mvec<T> p3ga2quadric2qc3ga(const p3ga2::Mvec<T> & mvQuadric){
    std::vector<T> quadCoeff = p3ga2::ga2quadric<T>(mvQuadric);
        return qc3ga::dualQuadric2ga<T>(quadCoeff);
    }

    /// \brief transform a qc3ga quadric to a p3ga2 one
    template<typename T>
    p3ga2::Mvec<T> qc3gaquadric2p3ga2(const qc3ga::Mvec<T> & mvQuadric){
    std::vector<T> quadCoeff = qc3ga::ga2DualQuadric<T>(mvQuadric);
        return p3ga2::quadric2ga<T>(quadCoeff);
    }

    /// \brief transform a p3ga2 quadric to a qc3ga one
    template<typename T>
    qc3ga::Mvec<T> p3ga2quadric2c3ga2(const p3ga2::Mvec<T> & mvQuadric){
    std::vector<T> quadCoeff = p3ga2::ga2quadric<T>(mvQuadric);
        return c3ga2::quadric2ga<T>(quadCoeff);
    }

    /// \brief transform a qc3ga quadric to a p3ga2 one
    template<typename T>
    p3ga2::Mvec<T> c3ga2quadric2p3ga2(const c3ga2::Mvec<T> & mvQuadric){
    std::vector<T> quadCoeff = c3ga2::ga2quadric<T>(mvQuadric);
        return p3ga2::quadric2ga<T>(quadCoeff);
    }




#endif // projection_inclusion_guard
