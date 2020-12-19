Garamon of project_namespace
=================

Garamon (Geometric Algebra Recursive and Adaptative Monster) is C++ library dedicated for project_namespace Geometric Algebra.

## Features
    * MIT Licence allows free use in all software (including GPL and commercial)
    * multi-platform (Windows, Linux, Unix, Mac)


Install
=======

## Dependencies
    * 'Eigen 3.3.4'  or more [Eigen](http://eigen.tuxfamily.org)
    * 'Cmake 3.10' or more
    * (optinal) 'Doxygen' if you want a documentation page

## Compiler tested
    * gcc 5.4.0
    * clang 4
    * apple-clang 900.0.39.2
    * MinGW 7.2.0
    * MSVC 19.14.26430.0 (Visual Studio 15.7.3)

## Install for Linux-Mac from the terminal
    * Check the dependencies
    * From Garamon root directory
        * 'mkdir build'
        * 'cd build'
        * 'cmake ..'
        * check that the cmake output has no errors
        * 'make'
        * the binary library is on the 'build' directory
        * (optional) 'sudo make install'
        * (optional) 'make html' to generate Doxygen documentation that can be found in 'build/doc/doc-doxygen/html/index.html'
    * you can find a sample code dedicated to project_namespace in the 'sample directory'


## Install for Windows with MinGW, using Windows Power Shell
    * Check the dependencies
    * From Garamon root directory
        * 'mkdir build'
        * 'cd build'
        * 'cmake -G "MinGW Makefiles" ..'
        * check that the cmake output has no errors
        * 'make'
        * the binary library is on the 'build' directory
        * (optional) 'mingw32-make install'
    * you can find a sample code dedicated to project_namespace in the 'sample directory'


Usage
=====

## Use of project_namespace

    * refer to the sample code in 'sample/src' directory and to the documentation or to the [cheatSheet.txt](cheatSheet.txt)


Notes
=====

## Authors
    * Vincent Nozick (Universite Paris-Est Marne-la-Vallee, France)
    * Stephane Breuils (Universite Paris-Est Marne-la-Vallee, France)

## Contact
    * vincent.nozick (at) u-pem.fr

## Thanks
    * special thanks to Lars Tingelstad for the Python binding

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



