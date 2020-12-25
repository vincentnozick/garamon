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

## Compiler tested
    * gcc 5.4.0
    * clang 4
    * apple-clang 900.0.39.2
    * MinGW 7.2.0
    * MSVC 19.14.26430.0 (Visual Studio 15.7.3)

## Install for Linux-Mac from the terminal
    * Check the dependencies
    * From Garamon Generator root directory
        * 'mkdir build'
        * 'cd build'
        * 'cmake ..'
        * check that the cmake output has no errors
        * 'make'
        * the binary executable is on the 'build' directory

## Install for Windows with MinGW, using Windows Power Shell
    * Check the dependencies
    * From Garamon Generator root directory
        * 'mkdir build'
        * 'cd build'
        * 'cmake -G "MinGW Makefiles" ..'
        * check that the cmake output has no errors
		* 'mingw32-make'
        * the binary executable is on the 'build' directory

## Install for Windows with Visual Studio 15 2017 Win64, using Windows Power Shell or cmd
    * Check the dependencies
    * From Garamon Generator root directory
        * 'mkdir build'
        * 'cd build'
        * 'cmake -G "Visual Studio 15 2017 Win64" ..'
        * check that the cmake output has no errors
		* open the file garamon_generator.sln with Visual Studio
		* generate the project "ALL_BUILD" with Release configuration
        * the binary executable is on the 'Release' directory


Usage
=====

## Generate a library

    * define the algebra to generate: chose a configuration file (.conf) on the 'conf' directory or create your own.
    * run the binary executable (from the 'build' directory) with the configuration file as argument
      > ./garamon_generator file.conf
    * the generated library is located in 'build/output' directory
    * to install the generated library, see its README.md

## Run the Python binding sample for a specific algebra
Let us consider the considered algebra is CGA of R3 corresponding to the configuration file c3ga.conf. 

    * Check the dependencies
    * From Garamon Generator root directory
    	* 'mkdir build'
    	* 'cd build'
    	* 'cmake ..'
    	* 'make'
    	* './garamon_generator ../conf/c3ga.conf'
    	* 'cd output/garamon_c3ga/'
    	* 'python setup.py build'
    	* 'python setup.py install'
    	* 'cd sample'
    	* 'python sample.py'

Notes
=====

## Authors
    * Vincent Nozick (Universite Paris-Est Marne-la-Vallee, France)
    * Stephane Breuils (National Institute of Informatics, Japan)

## Contact
    * vincent.nozick (at) u-pem.fr

## Reference
If you use Garamon for research purpose, please cite the following paper:

	@Article{breuils_garamon_2019,
	author="Breuils, St{\'e}phane and Nozick, Vincent and Fuchs, Laurent",
	title="Garamon: A Geometric Algebra Library Generator",
	journal="Advances in Applied Clifford Algebras",
	year="2019",
	month="Jul",
	day="22",
	volume="29",
	number="4",
	pages="69",
	issn="1661-4909",
	doi="10.1007/s00006-019-0987-7",
	url="https://doi.org/10.1007/s00006-019-0987-7"
	}

