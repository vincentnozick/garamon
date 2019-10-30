// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// p3gaTools.hpp
// Authors: Vincent Nozick and Stephane Breuils
// Contact: vincent.nozick@u-pem.fr

// These tools are dedicated to the dual version of R(3,0,1), meaning that the vectors are planes, bivectors are lines and vectors are points.
// For more details :
//
// Charles Gunn, Projective geometric algebra: A new framework for doing euclidean geometry, 2019.
// https://arxiv.org/pdf/1901.05873.pdf
//
// Charles Gunn, Course notes, Geometric Algebra for Computer Graphics, 2019.
// https://www.researchgate.net/publication/334397782_Course_notes_Geometric_Algebra_for_Computer_Graphics_SIGGRAPH_2019


/// \file p3gaTool.hpp
/// \author Stephane Breuils, Vincent Nozick
/// \brief some useful functions when using Projective Geometric Algebra of R(2,0,1). Use this file if you generated the lib using the "p3ga.conf" configuration file.


// Anti-doublon
#ifndef GARAMON_P3GA_TOOLS_HPP__
#define GARAMON_P3GA_TOOLS_HPP__
#pragma once


// External includes
#include <vector>
#include <random>
#include <chrono>

// Internal includes
#include <p3ga/Mvec.hpp>


// namespace for this specific algebra
namespace p3ga {


    /// \brief return the norm of a finite points or regular lines.
    /// \param mv corresponding to a finite point or regular line.
    /// \return sqrt(a.a + b.b) for lines and w for points
    template<typename T>
    T norm(const p3ga::Mvec<T>& mv) {
        return sqrt( mv | (~mv) );
    }

    /// \brief return the ideal norm of an ideal points or ideal lines.
    /// \param mv corresponding to an ideal points or ideal lines.
    /// \return 'c' for lines and sqrt(x.x + y.y) for points.
    template<typename T>
    T normIdeal(const p3ga::Mvec<T>& mv) {
        return norm(!mv);
    }

    /// \brief build a point from 3 components
    /// \param x vector component related to e023
    /// \param y vector component related to e013
    /// \param z vector component related to e012
    /// \return a multivector corresponding to a point p = e123 - x e023 + y e013 -z e012.
    template<typename T>
    p3ga::Mvec<T> point(const T &x, const T &y, const T &z) {

        p3ga::Mvec<T> mv;
        mv[p3ga::E023] = -x;
        mv[p3ga::E013] =  y;
        mv[p3ga::E012] = -z;
        mv[p3ga::E123] = 1.0;
        return mv;
    }

    // anonymous namespace for a random generator
    namespace {
        std::default_random_engine randomGenerator(0);
    }

    /// \brief reset the random generator seed
    /// \param seed used to initialize the random generator
    void randomSeed(const unsigned int seed){
        randomGenerator = std::default_random_engine(seed);
    }

    /// \brief reset the random generator seed using the time
    void randomSeed(){
        randomSeed(std::chrono::system_clock::now().time_since_epoch().count());
    }

    /// \brief build a random point with Euclidean coordinates ranging in [-1,1]
    /// \param seed optional random seem.
    /// \return a multivector corresponding to a point p = e123 - x e023 + y e013 -z e012.
    template<typename T>
    p3ga::Mvec<T> randomPoint(const unsigned seed=0){

        // uniform distribution over [-1,1]
        std::uniform_real_distribution<T> uniformRealDistribution(-1.0,1.0);

        // build the point
        return point<T>(uniformRealDistribution(randomGenerator), uniformRealDistribution(randomGenerator), uniformRealDistribution(randomGenerator));
    }

    /// \brief build a vector (ideal point / point at infinity) from 3 components
    /// \param x vector component related to e023
    /// \param y vector component related to e013
    /// \param z vector component related to e012
    /// \return a multivector corresponding to a point p = - x e02 + y e01
    template<typename T>
    p3ga::Mvec<T> vector(const T &x, const T &y, const T &z) {

        p3ga::Mvec<T> mv;
        mv[p3ga::E023] = -x;
        mv[p3ga::E013] =  y;
        mv[p3ga::E012] = -z;
        return mv;
    }

    template<typename T>
    p3ga::Mvec<T> lineSupportVector(const p3ga::Mvec<T> &line) {
        p3ga::Mvec<double> vec = p3ga::vector<T>(-line[E23], -line[E13], -line[E12]);
        return vec/normIdeal(vec);
    }

    /// \brief computes the nearest point from the origine that lies on the line (just project the origine on the line).
    /// \param line : the line to be computed.
    /// \return a point or an ideal point, not normalized.
    template<typename T>
    p3ga::Mvec<T> lineOrigine(const p3ga::Mvec<T> &line) {
        return (p3ga::e123<T>() | line) * line;
    }

    /// \brief build a plane with equation ax + by +cz +d = 0
    /// \param a line parameter related to x
    /// \param b line parameter related to y
    /// \param c line parameter related to z
    /// \param d line parameter (constant part)
    /// \return multivector : a.e1 + b.e2 + c.e3 + d.e0
    template<typename T>
    p3ga::Mvec<T> plane(const T &a, const T &b, const T&c, const T&d) {

        p3ga::Mvec<T> mv;
        mv[p3ga::E1] = a;
        mv[p3ga::E2] = b;
        mv[p3ga::E3] = c;
        mv[p3ga::E0] = d;
        return mv;
    }

    /// \brief return the join (antiwedge) of two multivectors
    /// \param mv1 left hand side multivector
    /// \param mv2 right hand side multivector
    /// \return mv1 v mv2 = !( !mv1 ^ !mv2 )
    template<typename T>
    p3ga::Mvec<T> operator&(const p3ga::Mvec<T>& mv1, const p3ga::Mvec<T>& mv2) {
        return ( mv1.dual() ^ mv2.dual() ).dual();
    }

} // namespace

#endif // projection_inclusion_guard