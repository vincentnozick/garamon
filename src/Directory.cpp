// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// Directory.cpp
// This file is part of the Garamon Generator.
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program


#include "Directory.hpp"

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <system_error>

namespace fs = std::filesystem;


void makeDirectory(const std::string &dirName) {

    // fails if the directory already exists or if its parent directory is missing (not recursive)
    std::error_code error;
    if(!fs::create_directory(dirName, error))
        throw std::runtime_error("can not create directory: " + dirName);
}


bool directoryOrFileExists(const std::string &dirName) {
    std::error_code error; // (an inaccessible path does not exist)
    return fs::exists(dirName, error);
}


std::string readFile(const std::string &fileName) {

    std::ifstream myfile;
    myfile.open(fileName, std::ios::in);

    // check if the file is opened
    if(!myfile.is_open()){
        throw std::runtime_error("can not open file: " + fileName);
    }

    // copy the data to a string
    std::string data;
    data.assign( (std::istreambuf_iterator<char>(myfile)) , (std::istreambuf_iterator<char>()) );

    myfile.close();

    return data;
}


bool writeFile(const std::string &data, const std::string &fileName) {

    std::ofstream myfile(fileName);
    if(!myfile.is_open()){
        std::cerr << "error: can not create file: " << fileName << std::endl;
        return false;
    }

    myfile << data;
    myfile.close();

    return true;
}





void substitute(std::string &data, const std::string &pattern, const std::string &replaceBy) {

    // literal replacement of all the occurrences of 'pattern' (neither 'pattern' nor 'replaceBy' are interpreted)
    if(pattern.empty())
        return;

    std::string result;
    result.reserve(data.size());
    std::size_t start = 0;
    std::size_t pos;
    while((pos = data.find(pattern, start)) != std::string::npos) {
        result.append(data, start, pos - start);
        result += replaceBy;
        start = pos + pattern.size();
    }
    result.append(data, start, std::string::npos);
    data.swap(result);
}


std::string joinPath(const std::string &directory, const std::string &name) {

    if(directory.empty())
        return name;

    const char last = directory.back();
    if(last == '/' || last == '\\')
        return directory + name;

    return directory + "/" + name;
}
