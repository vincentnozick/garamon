// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// ConfigParser.cpp
// This file is part of the Garamon Generator.
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program


#include "ConfigParser.hpp"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <iterator>

ConfigParser::ConfigParser(const std::string &filename) {

    // openfile in read only
    std::ifstream myfile;
    myfile.open(filename, std::ios::in);

    // check if the file is opened
    if(!myfile.is_open()){
        std::cerr << "error: can not open file: " << filename << std::endl;
        exit(EXIT_FAILURE);
    }else std::cout << "   open " << filename << " ... ok" << std::endl;

    // copy the data to a string
    data.assign( (std::istreambuf_iterator<char>(myfile)) , (std::istreambuf_iterator<char>()) );

    // close file
    myfile.close();
}

ConfigParser::~ConfigParser() {}

bool ConfigParser::extract(std::string &extractedData , const std::string &keyword) const {

    // find the starting point of the required data
    std::size_t start = data.find("<" + keyword + ">");

    // check if the keyword was found
    if(start == std::string::npos)
        return false;

    // read the end of the line (keyword + \n)
    start += ("<" + keyword + ">").size() + 1;

    // find the ending point of the required data
    std::size_t end = data.find("</" + keyword + ">", start);

    // remove the "\n"
    end--;

    // copy the data in the output string
    extractedData = data.substr(start, end-start);

    return true;
}

std::vector<double> vectorFromString(std::string const& stringData){
    std::istringstream iss(stringData);

    return std::vector<double>{
           std::istream_iterator<double>(iss),
           std::istream_iterator<double>()
    };
}

bool ConfigParser::readMatrix(const std::string &keyword, Eigen::MatrixXd &mat) const {

    std::vector<std::vector<double>> vector;

    // extract the matrix data from the file
    std::string stringTmp;
    if(!extract(stringTmp,keyword))
        return false;

    // split the string into a vector of vector
    std::string delimiter = " ";
    std::string lineDelimiter = "\n";
    vector.clear();
    size_t pos = 0;

    // for each line
    while ((pos = stringTmp.find(lineDelimiter)) != std::string::npos) {

        // remove duplicated delimiter
        if(pos == 0){
            stringTmp.erase(0, lineDelimiter.length());
            continue;
        }

        // extract a line
        std::string line = stringTmp.substr(0, pos);

        // add the new detected element
        vector.push_back(vectorFromString(line));

        // remove the extracted element from the input string
        stringTmp.erase(0, pos + delimiter.length());
    }

    // insert the last element (not followed by a delimiter)
    if(stringTmp.size() != 0)
        vector.push_back(vectorFromString(stringTmp));


    // convert the vectors into a matrix
    mat = Eigen::MatrixXd(vector.size(), vector[0].size());
    for(unsigned int i=0; i<(unsigned int)mat.rows(); ++i) {

        // if the current line size is not consistent with the matrix size
        if(vector[i].size() != (unsigned int) mat.cols())
            return false;

        for(unsigned int j = 0; j < (unsigned int)mat.cols(); ++j) {
            mat(i,j) = vector[i][j];
        }
    }

    return true;
}

bool ConfigParser::readString(const std::string &keyword, std::string &outputString) const {
    return extract(outputString, keyword);
}

bool ConfigParser::readUInt(const std::string &keyword, unsigned int &val) const {
    std::string stringTmp;
    val = 0;
    if(!extract(stringTmp,keyword))
        return false;
    val = std::stoi(stringTmp);
    return true;
}

bool ConfigParser::readDouble(const std::string &keyword, double &val) const {
    std::string stringTmp;
    val = 0;
    if(!extract(stringTmp,keyword))
        return false;
    val = std::stod(stringTmp);
    return true;
}

bool ConfigParser::readBool(const std::string &keyword, bool &val) const {
    std::string stringTmp;
    if(!extract(stringTmp,keyword))
        return false;

    // convert to lower case
    std::transform(stringTmp.begin(), stringTmp.end(), stringTmp.begin(), ::tolower);

    if(stringTmp == "true"){
        val = true;
        return true;
    }

    if(stringTmp == "false"){
        val = false;
        return true;
    }

    return false;
}

bool ConfigParser::readStringList(const std::string &keyword, std::vector<std::string> &stringList) const {

    // extract the string list from the string file
    std::string stringTmp;
    if(!extract(stringTmp,keyword))
        return false;

    // split the string into a vector of string
    std::string delimiter = " ";
    stringList.clear();
    size_t pos = 0;
    while ((pos = stringTmp.find(delimiter)) != std::string::npos) {

        // remove duplicated delimiter
        if(pos == 0){
            stringTmp.erase(0, delimiter.length());
            continue;
        }

        // add the new detected element
        stringList.push_back(stringTmp.substr(0, pos));

        // remove the extracted element from the input string
        stringTmp.erase(0, pos + delimiter.length());
    }

    // insert the last element (not followed by a delimiter)
    if(stringTmp.size() != 0)
        stringList.push_back(stringTmp);

    return true;
}


