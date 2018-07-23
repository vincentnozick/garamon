// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// c3ga2Tools.hpp
// Authors: Vincent Nozick and Stephane Breuils 
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


	
    /// Let us now build the extraction operators
    

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
    /// \brief construct the quadric whose coefficients are {a,b,c,d,e,f,g,h,i,j}
    template<typename T>
    std::vector<T> ga2quadric(const std::vector<T>& quadric){
        c3ga2::Mvec<T> mv = quadric[0]*Tx2<T>() + quadric[1]*Ty2<T>() + quadric[2]*Tz2<T>() + quadric[3]*Txy<T>() + quadric[4]*Tzx<T>() + quadric[5]*Tyz<T>()  
			  + quadric[6]*Tx<T>()  + quadric[7]*Ty<T>()  + quadric[8]*Tz<T>()  + quadric[9]*T1<T>();
	
        return mv;
    }



    /// \brief construct a Darboux cyclide with coefficients {a,b,c,d,e,f,g,h,i,j,k,l,m,n,o}, see Equation (53) of the paper
    /// Note that the same kind of function can be used to construct a Dupin cyclide
    template<typename T>
    std::vector<T> ga2DarbouxCyclide(const std::vector<T>& darbouxCyclide){
        c3ga2::Mvec<T> mv = darbouxCyclide[0]*Tt4<T>() + darbouxCyclide[2]*Txt2<T>() + darbouxCyclide[3]*Tyt2<T>() + darbouxCyclide[4]*Tzt2<T>() + darbouxCyclide[5]*Tt2<T>() + darbouxCyclide[6]*Tx2<T>() 
			  + darbouxCyclide[7]*Ty2<T>() + darbouxCyclide[8]*Tz2<T>() + darbouxCyclide[9]*Txy<T>() + darbouxCyclide[10]*Tzx<T>() + darbouxCyclide[11]*Tyz<T>()  
			  + darbouxCyclide[12]*Tx<T>()  + darbouxCyclide[13]*Ty<T>()  + darbouxCyclide[14]*Tz<T>()  + darbouxCyclide[15]*T1<T>();
        return mv;
    }







} // namespace

#endif // projection_inclusion_guard
