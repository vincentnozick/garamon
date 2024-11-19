// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// main.hpp
// This file is part of the Garamon Generator.
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

/// \file main.hpp
/// \author Stephane Breuils, Vincent Nozick
/// \brief Core of the library generator.


#include <iostream>
#include <algorithm>

#include "MetaData.hpp"
#include "Directory.hpp"
#include "Utilities.hpp"
#include "ProductToString.hpp"


int main(int argc, char** argv) {

    // Default output directory
    std::string outputDirectory = "output/";

    const std::string executableDir = getExecutableDirectory(argv[0]);
    std::string templateDataDirectory = executableDir + "/../data/";

    // Check the program arguments
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "-o") == 0) {
            if (i + 1 < argc) {
                outputDirectory = argv[++i];  // Get the next argument as the output directory
            } else {
                std::cerr << "error: option '-o' requires a directory path" << std::endl;
                return EXIT_FAILURE;
            }
        }
        if (std::strcmp(argv[i], "-d") == 0) {
            if (i + 1 < argc) {
                templateDataDirectory = argv[++i];  // Get the next argument as the data directory
            } else {
                std::cerr << "error: option '-d' requires a directory path" << std::endl;
                return EXIT_FAILURE;
            }
        }
    }

    // Check if the required configuration file argument is provided
    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " [-o output_dir] [-d data_dir] file.conf" << std::endl;
        std::cerr << "where 'file.conf' is the file that defines your algebra." << std::endl;
        return EXIT_FAILURE;
    }

    // Configuration file
    const std::string configurationFilesDirectory = argv[argc - 1];  // The last argument should be the config file



    // Read the configuration data
    std::cout << "load meta data ..." << std::endl;
    // Assuming MetaData class exists as shown in your code snippet
    MetaData metaData(configurationFilesDirectory);
    metaData.display();

    // Define the project directory structure
    std::cout << "define the project arborescence ..." << std::endl;
    std::string projectDirectory = outputDirectory + "/garamon_" + metaData.namespaceName;
    std::string srcDirectoryMain = projectDirectory + "/src";
    std::string srcDirectory = projectDirectory + "/src/" + metaData.namespaceName;
    std::string docDirectory = projectDirectory + "/doc";
    std::string docImageDirectory = projectDirectory + "/doc/images";
    std::string docHowToDirectory = projectDirectory + "/doc/HOWTO";
    std::string sampleDirectory = projectDirectory + "/sample";
    std::string srcSampleDirectory = projectDirectory + "/sample/src";
    std::string moduleSampleDirectory = projectDirectory + "/sample/modules";

    // Check if the directories are not already existing
    if (directoryOrFileExists(projectDirectory)) {
        std::cout.flush();  // Flushing the stream to avoid asynchronous output issues
        std::cerr << "error: a data or directory '" << projectDirectory << "' already exists" << std::endl;
        return EXIT_FAILURE;
    }

    // generate the arborescence at the right place
    makeDirectory(projectDirectory);
    makeDirectory(srcDirectoryMain);
    makeDirectory(srcDirectory);
    makeDirectory(docDirectory);
    makeDirectory(docImageDirectory);
    makeDirectory(docHowToDirectory);
    makeDirectory(sampleDirectory);
    makeDirectory(srcSampleDirectory);
    makeDirectory(moduleSampleDirectory);

    // namespace in uppercase for the #ifndef guards (antidoublons) and CMakeList.txt
    std::string upperCaseNamespace = metaData.namespaceName;
    std::transform(upperCaseNamespace.begin(), upperCaseNamespace.end(),upperCaseNamespace.begin(), ::toupper);

    // start to generate the files
    std::cout << "generate files ..." << std::endl;

    // add the fine files
    std::string data;

    // add the project CMakeLists.txt
    std::cout << "  cmake related files ..." << std::endl;
    data = readFile(templateDataDirectory + "CMakeLists.txt");
    substitute(data,"cmake_project_name_original_case", metaData.namespaceName);
    writeFile(data, projectDirectory + "/CMakeLists.txt");


    // add the documentation CMakeLists.txt
    data = readFile(templateDataDirectory + "doc/CMakeLists.txt");
    substitute(data,"cmake_project_name_original_case", metaData.namespaceName);
    writeFile(data, docDirectory + "/CMakeLists.txt");


    // add the documentation Doxyfile-html.cmake
    data = readFile(templateDataDirectory + "doc/Doxyfile-html.cmake");
    substitute(data,"cmake_project_name_original_case", metaData.namespaceName);
    if(metaData.dimension > 10) substitute(data,"cmake_project_ignore_constant_h", "*Constants.hpp");
    else substitute(data,"cmake_project_ignore_constant_h", "");
    writeFile(data, docDirectory + "/Doxyfile-html.cmake");


    // move the image
    //copyBin(templateDataDirectory + "doc/images/garamon.png", docImageDirectory + "/garamon.png");


    // mv HOWTO files
    copyText(templateDataDirectory + "doc/HOWTO/HOWTO-build", docHowToDirectory + "/HOWTO-build");
    copyText(templateDataDirectory + "doc/HOWTO/HOWTO-doc",   docHowToDirectory + "/HOWTO-doc");
    copyText(templateDataDirectory + "doc/HOWTO/HOWTO-test",  docHowToDirectory + "/HOWTO-test");


    // Initialize product tools:
    std::cout << "  data structure configuration ..." << std::endl;
    ProductTools algebraConfig(metaData.dimension);

    // bring some information on the grade, binomial coefficients and the index of any homogeneous multivector in the global multivector structure
    std::vector<int> perGradeStartingIndex = {0};
    computePerGradeStartingIndex(metaData.dimension, perGradeStartingIndex, 0, 0);

    // Basis transformation matrices initialization
    std::cout << "  basis transformation matrices initialization ..." << std::endl;
    std::vector<unsigned int> sizeTransformationMatrices; // for each matrix, number of non-zero elements
    std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor> > transformationMatrices;
    std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor> > inverseTransformationMatrices;

    // BasisTransformations.hpp
    std::cout << "  BasisTransformations.hpp ..." << std::endl;
    data = readFile(templateDataDirectory + "BasisTransformations.hpp");
    substitute(data,"project_inclusion_guard", upperCaseNamespace + "_BASISTRANSFORMATIONS_HPP__");
    substitute(data,"project_namespace", metaData.namespaceName);
    substitute(data,"project_algebra_dimension", std::to_string(metaData.dimension));
    if(metaData.inputMetricDiagonal == false) {
        // pair : matrix direct / inverse. These arrays contains triplet ([#row #col #coef])
        std::pair<std::vector<double>, std::vector<double> > allTransformationMatrices = computeTransformationMatricesToVector(
                metaData.transformationMatrix, metaData.epsilon, sizeTransformationMatrices, transformationMatrices, inverseTransformationMatrices);
        // put each transformation matrix data in a string
        substitute(data,"project_basischange_direct_loading", loadAllDirectOrInverseMatrices(sizeTransformationMatrices, allTransformationMatrices.first, false));
        substitute(data,"project_basischange_inverse_loading", loadAllDirectOrInverseMatrices(sizeTransformationMatrices, allTransformationMatrices.second, true));
        // sparse matrix initialization
        substitute(data,"project_calltobasischange_direct_loading", callAllDirectOrInverseMatricesFunctions(sizeTransformationMatrices, false) );
        substitute(data,"project_calltobasischange_inverse_loading", callAllDirectOrInverseMatricesFunctions(sizeTransformationMatrices, true) );
    }else{
        // diagonal metric = do not require any basis transformations
        substitute(data,"project_basischange_direct_loading", "");
        substitute(data,"project_basischange_inverse_loading","");
        substitute(data,"project_calltobasischange_direct_loading","" );
        substitute(data,"project_calltobasischange_inverse_loading","" );
    }
    writeFile(data, srcDirectory + "/BasisTransformations.hpp");


    // Constants.hpp
    std::cout << "  Constants.hpp ..." << std::endl;
    data = readFile(templateDataDirectory + "Constants.hpp");
    substitute(data,"project_inclusion_guard", upperCaseNamespace + "_CONSTANTS_HPP__");
    substitute(data,"project_namespace", metaData.namespaceName);
    substitute(data,"project_algebra_dimension", std::to_string(metaData.dimension));
    substitute(data,"project_basis_vectors_string", basisVectorsToString(metaData));
    substitute(data,"project_metric", metricToString(metaData));
    substitute(data,"project_per_grade_starting_index", perGradeStartingIndexToString(perGradeStartingIndex));
    substitute(data,"project_basis_vector_index", multivectorComponentBuilder(metaData,constantsDefinition()));
    substitute(data,"project_array_binomial_coefficient", binomialCoefToString(metaData.dimension));
    substitute(data,"project_sign_reverse", reverseSignArrayToString(metaData.dimension));
    substitute(data,"project_array_xorIndexConversion", xorIndexToGradeAndHomogeneousIndexArraysToString(metaData.dimension, algebraConfig));
    // the following string will contain the components of the fast dual array
    std::string fastDualComponents ="";
    if(metaData.fullRankMetric == false){ // no dual = replace by right complement
        substitute(data,"project_dual_arrays_permutations_and_coefficients", fastRightComplementUtilities(metaData.dimension,algebraConfig,srcDirectory,fastDualComponents));
        substitute(data,"project_dual_arrays_recursive_coefficients", "");
        substitute(data,"project_pseudo_scalar_inverse", "");
    }else{

        // compute the scale used to compute the inverse pseudo-scalar
        double signedPseudoScalar = (getScaleInversePseudoScalar(metaData.metric) * pow(-1,(metaData.dimension*(metaData.dimension-1))/2));
        substitute(data,"project_pseudo_scalar_inverse", "    constexpr double pseudoScalarInverse = " + std::to_string(signedPseudoScalar) + "; /*!< compute the inverse of the pseudo scalar */\n");       

        // fast dual
        if(metaData.inputMetricDiagonal == false) substitute(data,"project_dual_arrays_permutations_and_coefficients", fastDualUtilitiesBasisChange(metaData.dimension, algebraConfig, transformationMatrices, inverseTransformationMatrices, metaData.diagonalMetric,signedPseudoScalar,srcDirectory,fastDualComponents));
        else{substitute(data,"project_dual_arrays_permutations_and_coefficients", fastDualUtilities(metaData.dimension, algebraConfig,metaData.diagonalMetric,signedPseudoScalar,srcDirectory,fastDualComponents));}

        // recursive dual products
        if(metaData.inputMetricDiagonal == false) substitute(data,"project_dual_arrays_recursive_coefficients", primalWedgeDualUtilitiesBasisChange(metaData.dimension, algebraConfig,transformationMatrices,inverseTransformationMatrices,metaData.diagonalMetric,signedPseudoScalar) );
        else{substitute(data,"project_dual_arrays_recursive_coefficients", primalWedgeDualUtilities(metaData.dimension, algebraConfig, metaData.diagonalMetric,signedPseudoScalar) );}
    }
    substitute(data,"project_diagonal_Metric", diagonalMetricToString(metaData));
    if(metaData.inputMetricDiagonal == false) substitute(data,"project_load_transformation_matrices", basisTransformMatricesLoad());
    else{substitute(data,"project_load_transformation_matrices", "");}
    substitute(data,"project_dim_plus_one", std::to_string((metaData.dimension+1)));
    writeFile(data, srcDirectory + "/Constants.hpp");


    // DualCoefficients.hpp: load the components of the fast dual array
    data = readFile(templateDataDirectory + "DualCoefficients.hpp");
    substitute(data,"project_inclusion_guard", upperCaseNamespace + "_DUALCOEFFICIENTS_HPP__");
    substitute(data,"project_namespace", metaData.namespaceName);
    //if(metaData.fullRankMetric == false){ // no dual
    //    substitute(data,"project_fill_dual_array","");
    //}else{
    //    substitute(data,"project_fill_dual_array", loadAllDualCoefficientsArray(perGradeStartingIndex, fastDualComponents));
    // }
    substitute(data,"project_fill_dual_array", loadAllDualCoefficientsArray(perGradeStartingIndex, fastDualComponents));
    writeFile(data, srcDirectory + "/DualCoefficients.hpp");


    // Licence
    copyText(templateDataDirectory + "LICENCE.txt", projectDirectory + "/LICENCE.txt");


    // readme
    data = readFile(templateDataDirectory + "README.md");
    substitute(data,"project_namespace", metaData.namespaceName);
    writeFile(data, projectDirectory + "/README.md");


    // cheatSheet
    data = readFile(templateDataDirectory + "cheatSheet.txt");
    substitute(data,"project_namespace", metaData.namespaceName);
    substitute(data,"project_first_vector_basis", metaData.basisVectorName[0]);
    substitute(data,"project_second_vector_basis", metaData.basisVectorName[1]);
    if(bin_coeff(metaData.dimension,metaData.dimension/2) > metaData.maxDimPrecomputedProducts)  // if some constants are not created
        substitute(data,"project_limitation_vector_basis", "except for grade with basis vector of dimension higher than " + std::to_string(metaData.maxDimPrecomputedProducts));
    else substitute(data,"project_limitation_vector_basis", "");
    if(metaData.fullRankMetric) substitute(data,"project_limitation_dual", ""); // if dual exists
    else substitute(data,"project_limitation_dual", "\n------------------ dual are not defined in " + metaData.namespaceName + " ------------------"); // if dual does not exist
    writeFile(data, projectDirectory + "/cheatSheet.txt");


    // Utility.hpp
    data = readFile(templateDataDirectory + "Utility.hpp");
    substitute(data,"project_inclusion_guard", upperCaseNamespace + "_UTILITY_HPP__");
    substitute(data,"project_namespace", metaData.namespaceName);
    writeFile(data, srcDirectory + "/Utility.hpp");


    // Mvec.hpp
    std::cout << "  Mvec ..." << std::endl;
    data = readFile(templateDataDirectory + "Mvec.hpp");
    substitute(data,"project_inclusion_guard", upperCaseNamespace + "_MULTI_VECTOR_HPP__");
    substitute(data,"project_namespace", metaData.namespaceName);
    substitute(data,"project_multivector_one_component", multivectorComponentBuilder(metaData,oneComponentMultivectorPrototype())); // i.e. Mvec a = 2 * b.e12()
    substitute(data,"project_static_multivector_one_component", multivectorComponentBuilder(metaData,staticOneComponentMultivectorPrototype())); // i.e. Mvec a = 2 * cga::e12()
    substitute(data,"project_algebra_dimension", std::to_string(metaData.dimension));
    if(metaData.fullRankMetric == true){
        // keep the dual and other functions
        substitute(data,"project_singular_metric_comment_begin", "");
        substitute(data,"project_singular_metric_comment_end", "");
    }else{
        // comment the dual and other functions
        substitute(data,"project_singular_metric_comment_begin", singularMetricCommentBegin());
        substitute(data,"project_singular_metric_comment_end", singularMetricCommentEnd());       
    }
    substitute(data,"project_pseudo_scalar", std::to_string((1<<metaData.dimension)-1));
    // hybridization : maybe ignore sone products
    if(bin_coeff(metaData.dimension,metaData.dimension/2) > metaData.maxDimPrecomputedProducts){
        // only for geometric product, if the algebra dimension is too high to store all the precomputed products, insert a test to specify to use precomputed of recursive functions.
        substitute(data,"project_select_recursive_geometric_product_template", recursiveGeometricProductCallToString(metaData.maxDimPrecomputedProducts,!metaData.inputMetricDiagonal));
    }
    else  substitute(data,"project_select_recursive_geometric_product_template","");
    writeFile(data, srcDirectory + "/Mvec.hpp");


    // Mvec.cpp
    data = readFile(templateDataDirectory + "Mvec.cpp");
    substitute(data,"project_namespace", metaData.namespaceName);
    writeFile(data, srcDirectory + "/Mvec.cpp");


    // Outer.hpp
    std::cout << "  outer product ..." << std::endl;
    data = readFile(templateDataDirectory + "Outer.hpp");
    substitute(data,"project_inclusion_guard", upperCaseNamespace + "_OUTER_PRODUCT_HPP__");
    substitute(data,"project_namespace", metaData.namespaceName);
    if(metaData.fullRankMetric == true){
        // let the source code as it is
        substitute(data,"project_singular_metric_comment_begin", "");
        substitute(data,"project_singular_metric_comment_end", "");
    }else{
        // comment recursive outer product between dual
        substitute(data,"project_singular_metric_comment_begin", singularMetricCommentBegin());
        substitute(data,"project_singular_metric_comment_end", singularMetricCommentEnd());       
    }
    writeFile(data, srcDirectory + "/Outer.hpp");


    // OuterExplicit.hpp
    data = readFile(templateDataDirectory + "OuterExplicit.hpp");
    substitute(data,"project_inclusion_guard", upperCaseNamespace + "_OUTER_PRODUCT_EXPLICIT_HPP__");
    substitute(data,"project_namespace", metaData.namespaceName);
    substitute(data,"project_explicit_outer_functions", generateOuterExplicit_cpp(metaData,algebraConfig));
    substitute(data,"project_explicit_outer_pointer_functions", generateOuterExplicitFunctionsPointer(metaData.dimension));
    writeFile(data, srcDirectory + "/OuterExplicit.hpp");


    // Inner.hpp
    std::cout << "  inner product ..." << std::endl;
    if(metaData.identityMetric)  data = readFile(templateDataDirectory + "InnerEuclidean.hpp");
    else data = readFile(templateDataDirectory + "Inner.hpp");
    substitute(data,"project_inclusion_guard", upperCaseNamespace + "_INNER_PRODUCT_HPP__");
    substitute(data,"project_namespace", metaData.namespaceName);
    writeFile(data, srcDirectory + "/Inner.hpp");


    // InnerExplicit.hpp
    data = readFile(templateDataDirectory + "InnerExplicit.hpp");
    substitute(data,"project_inclusion_guard", upperCaseNamespace + "_INNER_PRODUCT_EXPLICIT_HPP__");
    substitute(data,"project_namespace", metaData.namespaceName);
    if(metaData.inputMetricDiagonal == false) substitute(data,"project_explicit_inner_functions", generateInnerExplicitBasisChange_cpp(metaData, transformationMatrices, inverseTransformationMatrices,algebraConfig));
    else{substitute(data,"project_explicit_inner_functions", generateInnerExplicit_cpp(metaData,algebraConfig));}
    substitute(data,"project_explicit_inner_pointer_functions", generateInnerExplicitFunctionsPointer(metaData.dimension));
    writeFile(data, srcDirectory + "/InnerExplicit.hpp");


    // Geometric.hpp
    std::cout << "  geometric product ..." << std::endl;
    if(metaData.identityMetric) data = readFile(templateDataDirectory + "GeometricEuclidean.hpp");
    else data = readFile(templateDataDirectory + "Geometric.hpp");
    substitute(data,"project_inclusion_guard", upperCaseNamespace + "_GEOMETRIC_PRODUCT_HPP__");
    substitute(data,"project_namespace", metaData.namespaceName);
    writeFile(data, srcDirectory + "/Geometric.hpp");


    // GeometricExplicit.hpp
    data = readFile(templateDataDirectory + "GeometricExplicit.hpp");
    substitute(data,"project_inclusion_guard", upperCaseNamespace + "_GEOMETRIC_PRODUCT_EXPLICIT_HPP__");
    substitute(data,"project_namespace", metaData.namespaceName);
    if(metaData.inputMetricDiagonal == false) substitute(data,"project_explicit_geometric_functions", generateGeometricExplicitBasisChange_cpp(metaData, transformationMatrices, inverseTransformationMatrices, algebraConfig));
    else{substitute(data,"project_explicit_geometric_functions", generateGeometricExplicit_cpp(metaData, algebraConfig));}
    substitute(data,"project_explicit_geometric_pointer_functions", generateGeometricExplicitFunctionsPointer(metaData));
    writeFile(data, srcDirectory + "/GeometricExplicit.hpp");


    // add the sample CMakeLists.txt
    std::cout << "  sample ..." << std::endl;
    data = readFile(templateDataDirectory + "sample/CMakeLists.txt");
    substitute(data,"cmake_project_name_sample", metaData.namespaceName + "_sample");
    substitute(data,"cmake_project_name_upper_case", upperCaseNamespace);
    substitute(data,"cmake_project_name_original_case", metaData.namespaceName);
    substitute(data,"cmake_project_name_original_case_py", metaData.namespaceName + "_py");
    writeFile(data, sampleDirectory + "/CMakeLists.txt");


    // add the sample modules 
    copyText(templateDataDirectory + "sample/modules/FindEigen.cmake", moduleSampleDirectory + "/FindEigen.cmake");
    data = readFile(templateDataDirectory + "sample/modules/Find_myLib.cmake");
    substitute(data,"cmake_project_name_upper_case", upperCaseNamespace);
    substitute(data,"cmake_project_name_original_case", metaData.namespaceName);
    writeFile(data, moduleSampleDirectory + "/Find" + upperCaseNamespace + ".cmake");


    // add the sample main.cpp
    data = readFile(templateDataDirectory + "sample/src/main.cpp");
    substitute(data,"project_namespace", metaData.namespaceName);
    substitute(data,"project_first_vector_basis", metaData.basisVectorName[0]);
    substitute(data,"project_second_vector_basis", metaData.basisVectorName[1]);
    writeFile(data, srcSampleDirectory + "/main.cpp");

    // PythonBindings.cpp
    data = readFile(templateDataDirectory + "PythonBindings.cpp");
    substitute(data,"project_namespace", metaData.namespaceName);
    substitute(data,"project_namespace_py", metaData.namespaceName + "_py");
    substitute(data,"project_static_multivector_one_component_python", multivectorComponentBuilder(metaData,staticOneComponentMultivectorPrototypePython())); // i.e. Mvec a = 2 * cga::e12()
    if(metaData.fullRankMetric == true){
        // keep the dual and other functions
        substitute(data,"project_singular_metric_comment_begin", "");
        substitute(data,"project_singular_metric_comment_end", "");
    }else{
        // comment the dual and other functions
        substitute(data,"project_singular_metric_comment_begin", singularMetricCommentBegin());
        substitute(data,"project_singular_metric_comment_end", singularMetricCommentEnd());       
    }
    substitute(data,"project_basis_vector_index", multivectorComponentBuilder(metaData,
    "m.attr(\"Eproject_name_blade\") = project_xor_index_blade;\n"
    ));
    writeFile(data, srcDirectory + "/PythonBindings.cpp");

    // setup.py
    data = readFile(templateDataDirectory + "setup.py");
    substitute(data,"project_namespace_py", metaData.namespaceName + "_py");
    writeFile(data, projectDirectory + "/setup.py");


    // sample.py
    data = readFile(templateDataDirectory + "sample/sample.py");
    substitute(data,"project_namespace_py", metaData.namespaceName + "_py");
    substitute(data,"project_first_vector_basis", metaData.basisVectorName[0]);
    substitute(data,"project_second_vector_basis", metaData.basisVectorName[1]);
    writeFile(data, sampleDirectory + "/sample.py");



    return EXIT_SUCCESS;
}

