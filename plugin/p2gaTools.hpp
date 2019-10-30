// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// p2gaTools.hpp
// Authors: Vincent Nozick and Stephane Breuils
// Contact: vincent.nozick@u-pem.fr

// These tools are dedicated to the dual version of R(2,0,1), meaning that the 1-vectors are lines and bivectors are points.
// For more details :
//
// Charles Gunn, Projective geometric algebra: A new framework for doing euclidean geometry, 2019.
// https://arxiv.org/pdf/1901.05873.pdf
//
// Charles Gunn, Course notes, Geometric Algebra for Computer Graphics, 2019.
// https://www.researchgate.net/publication/334397782_Course_notes_Geometric_Algebra_for_Computer_Graphics_SIGGRAPH_2019


/// \file p2gaTool.hpp
/// \author Stephane Breuils, Vincent Nozick
/// \brief some useful functions when using Projective Geometric Algebra of R(2,0,1). Use this file if you generated the lib using the "p2ga.conf" configuration file.


// Anti-doublon
#ifndef GARAMON_P2GA_TOOLS_HPP__
#define GARAMON_P2GA_TOOLS_HPP__
#pragma once


// External includes
#include <vector>
#include <random>
#include <chrono>

// Internal includes
#include <p2ga/Mvec.hpp>


// namespace for this specific algebra
namespace p2ga {

    /// \brief build a point from 2 components
    /// \param x vector component related to e02
    /// \param y vector component related to e01
    /// \return a multivector corresponding to a point p = e12 - x e02 + y e01
    template<typename T>
    p2ga::Mvec<T> point(const T &x, const T &y) {

        p2ga::Mvec<T> mv;
        mv[p2ga::E02] = -x;
        mv[p2ga::E01] =  y;
        mv[p2ga::E12] =  1.0;
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
    /// \return a multivector corresponding to a point p = e12 + x e02 + y e01.
    template<typename T>
    p2ga::Mvec<T> randomPoint(){

        // uniform distribution over [-1,1]
        std::uniform_real_distribution<T> uniformRealDistribution(-1.0,1.0);

        // build the random point
        return point<T>(uniformRealDistribution(randomGenerator), uniformRealDistribution(randomGenerator));
    }

    /// \brief build a vector (ideal point / point at infinity) from 2 components
    /// \param x vector component related to e02
    /// \param y vector component related to e01
    /// \return a multivector corresponding to a point p = - x e02 + y e01
    template<typename T>
    p2ga::Mvec<T> vector(const T &x, const T &y) {

        p2ga::Mvec<T> mv;
        mv[p2ga::E02] = -x;
        mv[p2ga::E01] =  y;
        return mv;
    }

    /// \brief build a line with equation ax + by +c = 0
    /// \param a line parameter related to x
    /// \param b line parameter related to y
    /// \param c line parameter (constant part)
    /// \return multivector : a.e1 + b.e2 + c.e0
    template<typename T>
    p2ga::Mvec<T> line(const T &a, const T &b, const T&c) {

        p2ga::Mvec<T> mv;
        mv[p2ga::E1] = a;
        mv[p2ga::E2] = b;
        mv[p2ga::E0] = c;
        return mv;
    }

    /// \brief return the join (antiwedge) of two multivectors
    /// \param mv1 left hand side multivector
    /// \param mv2 right hand side multivector
    /// \return mv1 v mv2 = !( !mv1 ^ !mv2 )
    template<typename T>
    p2ga::Mvec<T> operator&(const p2ga::Mvec<T>& mv1, const p2ga::Mvec<T>& mv2) {
        return ( mv1.dual() ^ mv2.dual() ).dual();
    }

    /// \brief return the norm of a finite points or regular lines.
    /// \param mv corresponding to a finite point or regular line.
    /// \return sqrt(a.a + b.b) for lines and w for points
    template<typename T>
    T norm(const p2ga::Mvec<T>& mv) {
        return sqrt( mv | (~mv) );
    }

    /// \brief return the ideal norm of an ideal points or ideal lines.
    /// \param mv corresponding to an ideal points or ideal lines.
    /// \return 'c' for lines and sqrt(x.x + y.y) for points.
    template<typename T>
    T normIdeal(const p2ga::Mvec<T>& mv) {
        return norm(!mv);
    }

    /// \brief computes the length of a closed loop polygon.
    /// \param polygon as a set of points.
    /// \return the length of the polygon.
    template<typename T>
    T closedLoopLength(const std::vector<p2ga::Mvec<T> > &polygon) {

        // compute the norm of each point
        std::vector<T> verticesNorm;
        for(const auto & vertex : polygon)
            verticesNorm.push_back(p2ga::norm(vertex));

        // compute the polygon length
        T sum = T(0);

        for(unsigned int i=0; i<polygon.size()-1; ++i)
            sum += p2ga::norm<T>( polygon[i] & polygon[i+1]) / (verticesNorm[i] * verticesNorm[i+1]);

        return sum;
    }

    /// \brief computes the area of a closed loop polygon.
    /// \param polygon as a set of points.
    /// \return the area of the polygon.
    template<typename T>
    T closedLoopArea(const std::vector<p2ga::Mvec<T> > &polygon) {

        // compute the norm of each point
        std::vector<T> verticesNorm;
        for(const auto & vertex : polygon)
            verticesNorm.push_back(p2ga::norm(vertex));

        // compute the polygon area
        p2ga::Mvec<double> sum;

        for(unsigned int i=0; i<polygon.size()-1; ++i)
            sum += (polygon[i] & polygon[i+1]) / (verticesNorm[i] * verticesNorm[i+1]);

        return T(0.5) * p2ga::normIdeal<T>(sum);
    }


} // namespace

#endif // projection_inclusion_guard