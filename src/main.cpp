// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// main.cpp
// This file is part of the Garamon Generator.
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

/// \file main.cpp
/// \author Stephane Breuils, Vincent Nozick
/// \brief Command line interface of the library generator (see Generator.hpp).


#include <iostream>
#include <cstdlib>
#include <exception>
#include <string>

#include <CLI/CLI.hpp>
#include "Generator.hpp"


int main(int argc, char** argv){

    // get the program arguments
    CLI::App app{"GARAMON: Geometric Algebra Recursive and Adaptative MONster Generator"};
    argv = app.ensure_utf8(argv);

    app.get_formatter()->column_width(35);
    app.get_formatter()->enable_option_type_names(false);
    app.get_formatter()->long_option_alignment_ratio( 0.2f /* 1/5*/);

    std::string configurationFilesDirectory;
    app.add_option("-c,--config", configurationFilesDirectory, "Configuration file, e.g. \"garamon/conf/c3ga.conf\"")
        ->check(CLI::ExistingFile)
        ->required();

    std::string templateDataDirectory;
    app.add_option("-t,--template", templateDataDirectory, "Template directory, e.g. \"garamon/data/\"")
        ->check(CLI::ExistingDirectory)
        ->required();

    std::string outputDirectory;
    app.add_option("-o,--output", outputDirectory, "Output directory, e.g. \"garamon/build/output/\"")
        ->check(CLI::ExistingDirectory)
        ->required();

    CLI11_PARSE(app, argc, argv);

    try {
        generateLibrary(configurationFilesDirectory, templateDataDirectory, outputDirectory);
    } catch(const std::exception &e) {
        std::cout.flush();
        std::cerr << "error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
