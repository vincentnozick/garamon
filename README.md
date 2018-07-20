Garamon Generator
=================

Garamon (Geometric Algebra Recursive and Adaptative Monster) is a generator of C++ libraries dedicated to Geometric Algebra.
From a configuration file, the library generator generate the source code of the specified Algebra as well as a compatible sample code and some documentation.

## Features
    * MIT Licence allows free use in all software (including GPL and commercial)
    * multi-platform (Windows, Linux, Unix, Mac)


Install
=======

## Dependencies
    * 'Eigen 3.3.4'  or more [Eigen](http://eigen.tuxfamily.org)
    * 'Cmake 3.10' or more
    * (optional) 'Doxygen' if you want documentation pages

## Compilator tested
    * gcc 5.4.0
    * clang 4
    * apple-clang 900.0.39.2
    * MinGW 7.2.0

## Install for Linux-Mac from the terminal
    * Check the dependencies
    * From Garamon Generator root directory
        * 'mkdir build'
        * 'cd build'
        * 'cmake ..'
        * check that the cmake output has no errors
        * 'make'
        * the binary executable is on the 'build' directory
        * (optional and not required) 'sudo make install'

## Install for Windows with MinGW, using Windows Power Shell
    * Check the dependencies
    * From Garamon Generator root directory
        * 'mkdir build'
        * 'cd build'
        * check that the cmake output has no errors
		* cmake -G "MinGW Makefiles" ..
		* mingw32-make
        * the binary executable is on the 'build' directory
        * (optional and not required) 'mingw32-make install'


Usage
=====

## Generate a library

    * define the algebra to generate: chose a configuration file (.conf) on the 'conf' directory or create your own.
    * run the binary executable (from the 'build' directory) with the configuration file as argument
      > ./garamon_generator file.conf
    * the generated library is located in 'build/output' directory


Notes
=====

## Authors
    * Vincent Nozick (Universite Paris-Est Marne-la-Vallee, France)
    * Stephane Breuils (Universite Paris-Est Marne-la-Vallee, France)

## Contact
    * vincent.nozick (at) u-pem.fr
