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
#include <cctype>
#include <cstddef>
#include <exception>
#include <fstream>
#include <ios>
#include <iostream>
#include <iterator>
#include <limits>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Core>

ConfigParser::ConfigParser(const std::string &filename) {

    // openfile in read only
    std::ifstream myfile;
    myfile.open(filename, std::ios::in);

    // check if the file is opened
    if(!myfile.is_open())
        throw std::runtime_error("can not open file: " + filename);
    std::cout << "   open " << filename << " ... ok" << std::endl;

    // copy the data to a string
    data.assign( (std::istreambuf_iterator<char>(myfile)) , (std::istreambuf_iterator<char>()) );

    // close file
    myfile.close();

    // files with Windows line endings: '\r' is just noise for this parser
    data.erase(std::remove(data.begin(), data.end(), '\r'), data.end());
}

ConfigParser::ConfigParser() {}

ConfigParser ConfigParser::fromString(const std::string &content) {
    ConfigParser parser;
    parser.data = content;
    parser.data.erase(std::remove(parser.data.begin(), parser.data.end(), '\r'), parser.data.end());
    return parser;
}

ConfigParser::~ConfigParser() {}

bool ConfigParser::extract(std::string &extractedData , const std::string &keyword) const {

    // find the starting point of the required data
    const std::string openingTag = "<" + keyword + ">";
    std::size_t start = data.find(openingTag);

    // check if the keyword was found
    if(start == std::string::npos)
        return false;
    start += openingTag.size();

    // find the ending point of the required data
    std::size_t end = data.find("</" + keyword + ">", start);
    if(end == std::string::npos)
        return false;

    // copy the data in the output string, without the surrounding blanks (the tags are usually on their own line)
    extractedData = data.substr(start, end-start);
    const char* blanks = " \t\n\r";
    std::size_t first = extractedData.find_first_not_of(blanks);
    if(first == std::string::npos) {
        extractedData.clear();
        return true;
    }
    std::size_t last = extractedData.find_last_not_of(blanks);
    extractedData = extractedData.substr(first, last - first + 1);

    return true;
}

// parse a line of numbers, return false if the line contains something else than numbers
static bool vectorFromString(std::string const& stringData, std::vector<double> &values){
    std::istringstream iss(stringData);
    values.clear();
    double value;
    while(iss >> value)
        values.push_back(value);
    return iss.eof();
}

bool ConfigParser::readMatrix(const std::string &keyword, Eigen::MatrixXd &mat) const {

    std::vector<std::vector<double>> vector;

    // extract the matrix data from the file
    std::string stringTmp;
    if(!extract(stringTmp,keyword))
        return false;

    // split the string into lines of numbers, ignoring the empty lines
    std::istringstream lines(stringTmp);
    std::string line;
    while(std::getline(lines, line)) {
        std::vector<double> row;
        if(!vectorFromString(line, row))
            return false;
        if(!row.empty())
            vector.push_back(row);
    }

    // empty matrix
    if(vector.empty())
        return false;

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

    // only digits are accepted (no sign, no decimal point, no trailing characters)
    if(stringTmp.empty() || !std::all_of(stringTmp.begin(), stringTmp.end(), [](unsigned char c){return std::isdigit(c) != 0;}))
        return false;

    try {
        unsigned long value = std::stoul(stringTmp);
        if(value > std::numeric_limits<unsigned int>::max())
            return false;
        val = (unsigned int) value;
    } catch(const std::exception&) {
        return false;
    }
    return true;
}

bool ConfigParser::readDouble(const std::string &keyword, double &val) const {
    std::string stringTmp;
    val = 0;
    if(!extract(stringTmp,keyword))
        return false;

    try {
        std::size_t parsed = 0;
        double value = std::stod(stringTmp, &parsed);
        if(parsed != stringTmp.size())
            return false;
        val = value;
    } catch(const std::exception&) {
        return false;
    }
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

    // split the string into a vector of string (any blank is a delimiter)
    stringList.clear();
    std::istringstream iss(stringTmp);
    std::string word;
    while(iss >> word)
        stringList.push_back(word);

    return true;
}
