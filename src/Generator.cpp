// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// Generator.cpp
// This file is part of the Garamon Generator.
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

/// \file Generator.cpp
/// \author Stephane Breuils, Vincent Nozick
/// \brief Core of the library generator.


#include <algorithm>
#include <cctype>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "Generator.hpp"
#include "MetaData.hpp"
#include "MetricTools.hpp"
#include "Directory.hpp"
#include "Utilities.hpp"
#include "ProductTools.hpp"
#include "ProductToString.hpp"


namespace {

    // write a generated file, throws if the file can not be written
    void writeGeneratedFile(const std::string &data, const std::string &fileName) {
        if(!writeFile(data, fileName))
            throw std::runtime_error("can not create file: " + fileName);
    }

    // copy a template file without modification, throws if the file can not be copied
    void copyTemplateFile(const std::string &src, const std::string &dest) {
        writeGeneratedFile(readFile(src), dest);
    }

}


std::string generateLibrary(const std::string &configurationFile,
                            const std::string &templateDataDirectory,
                            const std::string &outputDirectory){

    // read the configuration data
    std::cout << "load meta data ..." << std::endl;
    MetaData metaData(configurationFile);
    metaData.display();

    // template files
    auto templateFile = [&templateDataDirectory](const std::string &name){ return joinPath(templateDataDirectory, name); };

    // define the arborescence
    std::cout << "define the project arborescence ..." << std::endl;
    std::string projectDirectory      = joinPath(outputDirectory, "garamon_" + metaData.namespaceName);
    std::string srcDirectoryMain      = projectDirectory + "/src";
    std::string srcDirectory          = projectDirectory + "/src/" + metaData.namespaceName;
    std::string cmakeDirectory        = projectDirectory + "/cmake";
    std::string docDirectory          = projectDirectory + "/doc";
    std::string docImageDirectory     = projectDirectory + "/doc/images";
    std::string docHowToDirectory     = projectDirectory + "/doc/HOWTO";
    std::string sampleDirectory       = projectDirectory + "/sample";
    std::string srcSampleDirectory    = projectDirectory + "/sample/src";

    // check if the directories are not already existing
    if(directoryOrFileExists(projectDirectory))
        throw std::runtime_error("a data or directory '" + projectDirectory + "' already exists");

    // generate the arborescence at the right place
    makeDirectory(projectDirectory);
    makeDirectory(srcDirectoryMain);
    makeDirectory(srcDirectory);
    makeDirectory(cmakeDirectory);
    makeDirectory(docDirectory);
    makeDirectory(docImageDirectory);
    makeDirectory(docHowToDirectory);
    makeDirectory(sampleDirectory);
    makeDirectory(srcSampleDirectory);

    // namespace in uppercase for the #ifndef guards (antidoublons) and CMakeList.txt
    std::string upperCaseNamespace = metaData.namespaceName;
    std::transform(upperCaseNamespace.begin(), upperCaseNamespace.end(),upperCaseNamespace.begin(), ::toupper);

    // start to generate the files
    std::cout << "generate files ..." << std::endl;

    // add the fine files
    std::string data;

    // add the project CMakeLists.txt
    std::cout << "  cmake related files ..." << std::endl;
    data = readFile(templateFile("CMakeLists.txt"));
    substitute(data,"cmake_project_name_original_case", metaData.namespaceName);
    writeGeneratedFile(data, projectDirectory + "/CMakeLists.txt");


    // add the CMake package configuration file (find_package(<namespace> CONFIG))
    data = readFile(templateFile("cmake/Config.cmake.in"));
    substitute(data,"cmake_project_name_original_case", metaData.namespaceName);
    writeGeneratedFile(data, cmakeDirectory + "/" + metaData.namespaceName + "Config.cmake.in");


    // add the documentation CMakeLists.txt
    data = readFile(templateFile("doc/CMakeLists.txt"));
    substitute(data,"cmake_project_name_original_case", metaData.namespaceName);
    writeGeneratedFile(data, docDirectory + "/CMakeLists.txt");


    // add the documentation Doxyfile-html.cmake
    data = readFile(templateFile("doc/Doxyfile-html.cmake"));
    substitute(data,"cmake_project_name_original_case", metaData.namespaceName);
    if(metaData.dimension > 10) substitute(data,"cmake_project_ignore_constant_h", "*Constants.hpp");
    else substitute(data,"cmake_project_ignore_constant_h", "");
    writeGeneratedFile(data, docDirectory + "/Doxyfile-html.cmake");


    // mv HOWTO files
    copyTemplateFile(templateFile("doc/HOWTO/HOWTO-build"), docHowToDirectory + "/HOWTO-build");
    copyTemplateFile(templateFile("doc/HOWTO/HOWTO-doc"),   docHowToDirectory + "/HOWTO-doc");
    copyTemplateFile(templateFile("doc/HOWTO/HOWTO-test"),  docHowToDirectory + "/HOWTO-test");


    // Initialize product tools:
    std::cout << "  data structure configuration ..." << std::endl;
    ProductTools algebraConfig(metaData.dimension);

    // bring some information on the grade, binomial coefficients and the index of any homogeneous multivector in the global multivector structure
    std::vector<int> perGradeStartingIndex = {0};
    computePerGradeStartingIndex(metaData.dimension, perGradeStartingIndex, 0, 0);

    // Basis transformation matrices initialization
    std::cout << "  basis transformation matrices initialization ..." << std::endl;
    std::vector<unsigned int> sizeTransformationMatrices;        // for each matrix, number of non-zero elements
    std::vector<unsigned int> sizeInverseTransformationMatrices; // for each inverse matrix, number of non-zero elements
    std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor> > transformationMatrices;
    std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor> > inverseTransformationMatrices;

    // BasisTransformations.hpp
    std::cout << "  BasisTransformations.hpp ..." << std::endl;
    data = readFile(templateFile("BasisTransformations.hpp"));
    substitute(data,"project_inclusion_guard", upperCaseNamespace + "_BASISTRANSFORMATIONS_HPP__");
    substitute(data,"project_namespace", metaData.namespaceName);
    substitute(data,"project_algebra_dimension", std::to_string(metaData.dimension));
    if(metaData.inputMetricDiagonal == false) {
        // pair : matrix direct / inverse. These arrays contains triplet ([#row #col #coef])
        std::pair<std::vector<double>, std::vector<double> > allTransformationMatrices = computeTransformationMatricesToVector(
                metaData.transformationMatrix, metaData.epsilon, sizeTransformationMatrices, sizeInverseTransformationMatrices, transformationMatrices, inverseTransformationMatrices);
        // put each transformation matrix data in a string
        substitute(data,"project_basischange_direct_loading", loadAllDirectOrInverseMatrices(sizeTransformationMatrices, allTransformationMatrices.first, false));
        substitute(data,"project_basischange_inverse_loading", loadAllDirectOrInverseMatrices(sizeInverseTransformationMatrices, allTransformationMatrices.second, true));
        // sparse matrix initialization
        substitute(data,"project_calltobasischange_direct_loading", callAllDirectOrInverseMatricesFunctions(sizeTransformationMatrices, false) );
        substitute(data,"project_calltobasischange_inverse_loading", callAllDirectOrInverseMatricesFunctions(sizeInverseTransformationMatrices, true) );
    }else{
        // diagonal metric = do not require any basis transformations
        substitute(data,"project_basischange_direct_loading", "");
        substitute(data,"project_basischange_inverse_loading","");
        substitute(data,"project_calltobasischange_direct_loading","" );
        substitute(data,"project_calltobasischange_inverse_loading","" );
    }
    writeGeneratedFile(data, srcDirectory + "/BasisTransformations.hpp");


    // Constants.hpp
    std::cout << "  Constants.hpp ..." << std::endl;
    data = readFile(templateFile("Constants.hpp"));
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
        substitute(data,"project_pseudo_scalar_inverse", "");
    }else{

        // compute the scale used to compute the inverse pseudo-scalar
        double signedPseudoScalar = (getScaleInversePseudoScalar(metaData.metric) * std::pow(-1,(metaData.dimension*(metaData.dimension-1))/2));
        substitute(data,"project_pseudo_scalar_inverse", "    constexpr double pseudoScalarInverse = " + doubleToString(signedPseudoScalar) + "; /*!< compute the inverse of the pseudo scalar */\n");

        // fast dual
        if(metaData.inputMetricDiagonal) substitute(data,"project_dual_arrays_permutations_and_coefficients", fastDualUtilities(metaData.dimension, algebraConfig,metaData.diagonalMetric,signedPseudoScalar,srcDirectory,fastDualComponents));
        else if(metaData.inputMetricPermutationOfDiagonal) substitute(data,"project_dual_arrays_permutations_and_coefficients", fastDualUtilitiesBasisChange(metaData.dimension, algebraConfig, transformationMatrices, inverseTransformationMatrices, metaData.diagonalMetric,signedPseudoScalar,srcDirectory,fastDualComponents));
        else{
            // the dual of a basis blade is not a single basis blade: no fast dual (see fastDualAvailable), the arrays are only
            // generated with a consistent size (the right complement arrays) to keep the library code unchanged
            substitute(data,"project_dual_arrays_permutations_and_coefficients", fastRightComplementUtilities(metaData.dimension,algebraConfig,srcDirectory,fastDualComponents));
        }
    }
    // the fast dual (permutation and scaling) is exact for the right complement (degenerate metric), and for the dual if the metric
    // is diagonal or a permutation of a diagonal matrix (then the dual of a basis blade is a single basis blade)
    const bool fastDualAvailable = !metaData.fullRankMetric || metaData.inputMetricDiagonal || metaData.inputMetricPermutationOfDiagonal;
    substitute(data,"project_fast_dual_available", fastDualAvailable ? "true" : "false");
    substitute(data,"project_diagonal_Metric", diagonalMetricToString(metaData));
    if(metaData.inputMetricDiagonal == false) substitute(data,"project_load_transformation_matrices", basisTransformMatricesLoad());
    else{substitute(data,"project_load_transformation_matrices", "");}
    substitute(data,"project_dim_plus_one", std::to_string((metaData.dimension+1)));
    writeGeneratedFile(data, srcDirectory + "/Constants.hpp");


    // DualCoefficients.hpp: load the components of the fast dual array
    data = readFile(templateFile("DualCoefficients.hpp"));
    substitute(data,"project_inclusion_guard", upperCaseNamespace + "_DUALCOEFFICIENTS_HPP__");
    substitute(data,"project_namespace", metaData.namespaceName);
    substitute(data,"project_fill_dual_array", loadAllDualCoefficientsArray(perGradeStartingIndex, fastDualComponents));
    writeGeneratedFile(data, srcDirectory + "/DualCoefficients.hpp");


    // Licence
    copyTemplateFile(templateFile("LICENCE.txt"), projectDirectory + "/LICENCE.txt");


    // readme
    data = readFile(templateFile("README.md"));
    substitute(data,"project_namespace", metaData.namespaceName);
    writeGeneratedFile(data, projectDirectory + "/README.md");


    // cheatSheet
    data = readFile(templateFile("cheatSheet.txt"));
    substitute(data,"project_namespace", metaData.namespaceName);
    substitute(data,"project_first_vector_basis", metaData.basisVectorName[0]);
    substitute(data,"project_second_vector_basis", metaData.basisVectorName[metaData.dimension > 1 ? 1 : 0]);
    if(bin_coeff(metaData.dimension,metaData.dimension/2) > metaData.maxDimBasisAccessor)  // if some constants are not created
        substitute(data,"project_limitation_vector_basis", "except for grade with basis vector of dimension higher than " + std::to_string(metaData.maxDimBasisAccessor));
    else substitute(data,"project_limitation_vector_basis", "");
    if(metaData.fullRankMetric) substitute(data,"project_limitation_dual", ""); // if dual exists
    else substitute(data,"project_limitation_dual", "\n------------------ dual are not defined in " + metaData.namespaceName + " (degenerate metric): dual() computes the right complement ------------------"); // if dual does not exist
    writeGeneratedFile(data, projectDirectory + "/cheatSheet.txt");


    // Utility.hpp
    data = readFile(templateFile("Utility.hpp"));
    substitute(data,"project_inclusion_guard", upperCaseNamespace + "_UTILITY_HPP__");
    substitute(data,"project_namespace", metaData.namespaceName);
    writeGeneratedFile(data, srcDirectory + "/Utility.hpp");


    // Mvec.hpp
    std::cout << "  Mvec ..." << std::endl;
    data = readFile(templateFile("Mvec.hpp"));
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
    writeGeneratedFile(data, srcDirectory + "/Mvec.hpp");


    // Outer.hpp
    std::cout << "  outer product ..." << std::endl;
    data = readFile(templateFile("Outer.hpp"));
    substitute(data,"project_inclusion_guard", upperCaseNamespace + "_OUTER_PRODUCT_HPP__");
    substitute(data,"project_namespace", metaData.namespaceName);
    writeGeneratedFile(data, srcDirectory + "/Outer.hpp");


    // OuterExplicit.hpp
    data = readFile(templateFile("OuterExplicit.hpp"));
    substitute(data,"project_inclusion_guard", upperCaseNamespace + "_OUTER_PRODUCT_EXPLICIT_HPP__");
    substitute(data,"project_namespace", metaData.namespaceName);
    substitute(data,"project_explicit_outer_functions", generateOuterExplicit_cpp(metaData,algebraConfig));
    substitute(data,"project_explicit_outer_pointer_functions", generateOuterExplicitFunctionsPointer(metaData.dimension));
    writeGeneratedFile(data, srcDirectory + "/OuterExplicit.hpp");


    // Inner.hpp
    std::cout << "  inner product ..." << std::endl;
    if(metaData.identityMetric)  data = readFile(templateFile("InnerEuclidean.hpp"));
    else data = readFile(templateFile("Inner.hpp"));
    substitute(data,"project_inclusion_guard", upperCaseNamespace + "_INNER_PRODUCT_HPP__");
    substitute(data,"project_namespace", metaData.namespaceName);
    writeGeneratedFile(data, srcDirectory + "/Inner.hpp");


    // InnerExplicit.hpp
    data = readFile(templateFile("InnerExplicit.hpp"));
    substitute(data,"project_inclusion_guard", upperCaseNamespace + "_INNER_PRODUCT_EXPLICIT_HPP__");
    substitute(data,"project_namespace", metaData.namespaceName);
    if(metaData.inputMetricDiagonal == false) substitute(data,"project_explicit_inner_functions", generateInnerExplicitBasisChange_cpp(metaData, transformationMatrices, inverseTransformationMatrices,algebraConfig));
    else{substitute(data,"project_explicit_inner_functions", generateInnerExplicit_cpp(metaData,algebraConfig));}
    substitute(data,"project_explicit_inner_pointer_functions", generateInnerExplicitFunctionsPointer(metaData.dimension));
    writeGeneratedFile(data, srcDirectory + "/InnerExplicit.hpp");


    // Geometric.hpp
    std::cout << "  geometric product ..." << std::endl;
    if(metaData.identityMetric) data = readFile(templateFile("GeometricEuclidean.hpp"));
    else data = readFile(templateFile("Geometric.hpp"));
    substitute(data,"project_inclusion_guard", upperCaseNamespace + "_GEOMETRIC_PRODUCT_HPP__");
    substitute(data,"project_namespace", metaData.namespaceName);
    writeGeneratedFile(data, srcDirectory + "/Geometric.hpp");


    // GeometricExplicit.hpp
    data = readFile(templateFile("GeometricExplicit.hpp"));
    substitute(data,"project_inclusion_guard", upperCaseNamespace + "_GEOMETRIC_PRODUCT_EXPLICIT_HPP__");
    substitute(data,"project_namespace", metaData.namespaceName);
    if(metaData.inputMetricDiagonal == false) substitute(data,"project_explicit_geometric_functions", generateGeometricExplicitBasisChange_cpp(metaData, transformationMatrices, inverseTransformationMatrices, algebraConfig));
    else{substitute(data,"project_explicit_geometric_functions", generateGeometricExplicit_cpp(metaData, algebraConfig));}
    substitute(data,"project_explicit_geometric_pointer_functions", generateGeometricExplicitFunctionsPointer(metaData));
    writeGeneratedFile(data, srcDirectory + "/GeometricExplicit.hpp");


    // add the sample CMakeLists.txt
    std::cout << "  sample ..." << std::endl;
    data = readFile(templateFile("sample/CMakeLists.txt"));
    substitute(data,"cmake_project_name_sample", metaData.namespaceName + "_sample");
    substitute(data,"cmake_project_name_original_case", metaData.namespaceName);
    writeGeneratedFile(data, sampleDirectory + "/CMakeLists.txt");


    // add the sample main.cpp
    data = readFile(templateFile("sample/src/main.cpp"));
    substitute(data,"project_namespace", metaData.namespaceName);
    substitute(data,"project_first_vector_basis", metaData.basisVectorName[0]);
    substitute(data,"project_second_vector_basis", metaData.basisVectorName[metaData.dimension > 1 ? 1 : 0]);
    writeGeneratedFile(data, srcSampleDirectory + "/main.cpp");

    // PythonBindings.cpp
    data = readFile(templateFile("PythonBindings.cpp"));
    substitute(data,"project_namespace_py", metaData.namespaceName + "_py");
    substitute(data,"project_namespace", metaData.namespaceName);
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
    writeGeneratedFile(data, srcDirectory + "/PythonBindings.cpp");

    // sample.py
    data = readFile(templateFile("sample/sample.py"));
    substitute(data,"project_namespace_py", metaData.namespaceName + "_py");
    substitute(data,"project_first_vector_basis", metaData.basisVectorName[0]);
    substitute(data,"project_second_vector_basis", metaData.basisVectorName[metaData.dimension > 1 ? 1 : 0]);
    writeGeneratedFile(data, sampleDirectory + "/sample.py");

    return projectDirectory;
}
