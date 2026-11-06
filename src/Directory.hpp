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

void makeDirectory(const std::string &dirName);

bool directoryOrFileExists(const std::string &dirName);

std::string readFile(const std::string& fileName);

bool writeFile(const std::string& data, const std::string& fileName);

// replace all the occurrences of the literal string 'pattern' in 'data' by the literal string 'replaceBy'
void substitute(std::string &data, const std::string &pattern, const std::string &replaceBy);

// concatenate a directory and a file/directory name, adding a '/' separator if needed (i.e. "dir" + "file" => "dir/file", "dir/" + "file" => "dir/file")
std::string joinPath(const std::string &directory, const std::string &name);


#endif  // GARAGEN_DIRECTORY_HPP__