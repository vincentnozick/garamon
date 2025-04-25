// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// Directory.hpp
// This file is part of the Garamon Generator.
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

/// \file Directory.hpp
/// \author Stephane Breuils, Vincent Nozick
/// \brief Basic tools to manipulate files and directories.


#ifndef GARAGEN_DIRECTORY_HPP__
#define GARAGEN_DIRECTORY_HPP__

#include <string>
#include <vector>

void makeDirectory(const std::string &dirName);

bool directoryExists(const std::string &dirName);

bool directoryOrFileExists(const std::string &dirName);

bool directoryOrFileExists_ifstream(const std::string& name);

std::string readFile(const std::string& fileName);

bool writeFile(const std::string& data, const std::string& fileName);

void substitute(std::string &data, const std::string &pattern, const std::string &replaceBy);

bool copyBin(const std::string &src, const std::string &dest);

bool copyText(const std::string &src, const std::string &dest);

std::string getExecutableDirectory(const char* executable);

#endif  // GARAGEN_DIRECTORY_HPP__