Garamon Generator [![Discord chat](https://img.shields.io/discord/1009227277491519618.svg?logo=discord&colorB=7289DA&style=flat-square)](https://discord.gg/dYe2bPAWEQ)
=================

Garamon (Geometric Algebra Recursive and Adaptative Monster) is a generator of C++ libraries dedicated to Geometric Algebra.
From a configuration file, the library generator generate the source code of the specified Algebra as well as a compatible sample code and some documentation.

For help or discussion related to Garamon, you can join us on [Discord](https://discord.gg/dYe2bPAWEQ).

## Features
    * MIT Licence allows free use in all software (including GPL and commercial)
    * multi-platform (Windows, Linux, Unix, Mac)


Install
=======

## Dependencies
    * 'Eigen 3.3.4'  or more [Eigen](http://eigen.tuxfamily.org)
    * 'CLI11 2.7.2'  or more [CLI11](https://github.com/CLIUtils/CLI11)
    * 'CMake 3.24' or more
    * 'Catch2 3' or more [Catch2](https://github.com/catchorg/Catch2) (only for the unit tests, downloaded if not found)
    * 'pybind11' (only for building the Python bindings of a generated library, install with 'pip install pybind11'; requires CMake 3.18 or more)

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

```sh
~/garamon/build$ ./garamon_generator --help
GARAMON: Geometric Algebra Recursive and Adaptative MONster Generator

./garamon_generator [OPTIONS]

OPTIONS:
  -h,  --help                      Print this help message and exit
  -c,  --config REQUIRED           Configuration file, e.g. "garamon/conf/c3ga.conf"
  -t,  --template REQUIRED         Template directory, e.g. "garamon/data/"
  -o,  --output REQUIRED           Output directory, e.g. "garamon/build/output/"

~/garamon/build$
```

## Generate a library

    * define the algebra to generate: chose a configuration file (.conf) on the 'conf' directory or create your own.
    * run the binary executable (from the 'build' directory) with the configuration file as argument
      > ./garamon_generator -c ../conf/ga.conf -t ../data/ -o ./output/
    * the generated library is located in 'build/output' directory
    * to install the generated library, see its README.md

## Run the Python binding sample for a specific algebra
Let us consider the considered algebra is CGA of R3 corresponding to the configuration file c3ga.conf. 

    * Check the dependencies (in particular `pip install pybind11`, Python 3 development headers and CMake 3.18 or more)
    * From Garamon Generator root directory
    	* 'mkdir build'
    	* 'cd build'
    	* 'cmake ..'
    	* 'make'
        * './garamon_generator -c ../conf/c3ga.conf -t ../data/ -o ./output/'
    	* 'cd output/garamon_c3ga/'
    	* 'cmake -S . -B build -DBUILD_PYTHON=ON' (¹) (²)
    	* 'cmake --build build --config Release'
    	* 'cmake --build build --config Release --target python_sample' (runs 'sample/sample.py' with the just-built module)
    	* (optional) 'cmake --install build --config Release --component python' (²) to make the module importable from anywhere, then:
    		* 'cd sample'
    		* 'python sample.py'

(1) The `cmake -S ...` command requires the Eigen3 library. It will fail if it can't find it.
You might need to add the `-DCMAKE_PREFIX_PATH=/path/to/eigen3` parameter to point to where the library is.

(2) The module is installed into the site-packages of the Python interpreter found by CMake, which usually requires admin rights (`sudo`, or an elevated terminal on Windows).
Either install to your per-user site-packages instead, by adding `-DPYTHON_USER_INSTALL=ON` to the `cmake -S ...` command,
or activate a virtual environment (`python -m venv .venv`) before the `cmake -S ...` command: CMake then uses it and installs the module there.


Tests
=====

The unit tests use [Catch2](https://github.com/catchorg/Catch2) (v3) and CTest. They are built by default (CMake option `GARAMON_BUILD_TESTS`):

```sh
~/garamon$ cmake -S . -B build
~/garamon$ cmake --build build
~/garamon$ ctest --test-dir build --output-on-failure
```

* `generator_tests`: unit tests of the generator itself (configuration parser, metric analysis, product tables, code generation, ...)
* `<algebra>_tests`: the algebras `c3ga`, `c5ga`, `e3ga` and `p3ga`, as well as two test algebras with a non-normalized and a non-orthogonal metric (`tests/conf/`), are generated at build time; each generated library is compared with an independent reference implementation of the geometric algebra (`tests/algebras/reference/CliffordReference.hpp`)
* `<algebra>: sample`: the generated sample program is compiled and run.

To run a part of the tests, e.g. the c3ga dual tests: `build/tests/c3ga_tests "[dual]"`, or `ctest --test-dir build -R c3ga`.


Notes
=====

## Authors
    * Vincent Nozick (Universite Paris-Est Marne-la-Vallee, France)
    * Stephane Breuils (Universite Savoie Mont-Blanc, Japan)

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

