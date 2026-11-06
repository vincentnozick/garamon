// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// Generator.hpp
// This file is part of the Garamon Generator.
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

/// \file Generator.hpp
/// \author Stephane Breuils, Vincent Nozick
/// \brief Core of the library generator.

#ifndef GARAGEN_GENERATOR_HPP__
#define GARAGEN_GENERATOR_HPP__

#include <string>


/// \brief generate the library of the algebra defined in a configuration file.
/// \param configurationFile - configuration file of the algebra, e.g. "garamon/conf/c3ga.conf"
/// \param templateDataDirectory - directory of the template files, e.g. "garamon/data/"
/// \param outputDirectory - existing directory where the library directory "garamon_[namespace]" is created, e.g. "garamon/build/output/"
/// \return the directory of the generated library, i.e. "[outputDirectory]/garamon_[namespace]"
/// \throws std::runtime_error if the configuration is invalid, if a file can not be read or written, or if the library directory already exists
std::string generateLibrary(const std::string &configurationFile,
                            const std::string &templateDataDirectory,
                            const std::string &outputDirectory);


#endif // GARAGEN_GENERATOR_HPP__
