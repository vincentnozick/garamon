// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// ConfigParser.hpp
// This file is part of the Garamon Generator.
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

/// \file ConfigParser.hpp
/// \author Stephane Breuils, Vincent Nozick
/// \brief Very basic parser.

#ifndef GARAGEN_CONFIG_PARSER_HPP__
#define GARAGEN_CONFIG_PARSER_HPP__

#include <iostream>
#include <vector>
#include <string>

#include <Eigen/Core>


class ConfigParser {

protected:
    std::string data;

public:
    ConfigParser(const std::string &filename);
    ~ConfigParser();

    bool extract(std::string &extractedData , const std::string &keyword) const;
    bool readMatrix(const std::string &keyword, Eigen::MatrixXd &mat) const;
    bool readString(const std::string &keyword, std::string &outputString) const;
    bool readUInt(const std::string &keyword, unsigned int &val) const;
    bool readDouble(const std::string &keyword, double &val) const;
    bool readBool(const std::string &keyword, bool &val) const;
    bool readStringList(const std::string &keyword, std::vector<std::string> &stringList) const;
};


#endif // GARAGEN_CONFIG_PARSER_HPP__