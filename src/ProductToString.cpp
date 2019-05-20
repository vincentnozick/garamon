// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// ProductToString.cpp
// This file is part of the Garamon Generator.
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

#include "ProductToString.hpp"
#include "Utilities.hpp"  // for binomial
#include "Directory.hpp"  // for substitute



// generate the name of the vector basis defined in the metaDat, used for the outputs, i.e. something like {"0", "1", "2", "3", "i"}
std::string basisVectorsToString(const MetaData &metaData){

    std::string basisList = "{";
    bool first = true;
    for(auto &i: metaData.basisVectorName){
        if(!first)
            basisList += ", ";
        basisList += "\"" + i + "\"";
        first = false;
    }

    basisList += "}";
    return basisList;
}


// generate a string containing an array specifying the index of each first element of grade k in the full multivector.
// input: the index array to convert to string
// example (c3ga) : static constexpr unsigned int perGradeStartingIndex[6] = {0,1,6,16,26,31};
std::string perGradeStartingIndexToString(const std::vector<int>& perGradeStartingIndex){

    std::string indicesList = "perGradeStartingIndex["+std::to_string(perGradeStartingIndex.size())+"] = {";

    for(size_t i=0;i<perGradeStartingIndex.size();++i){
        indicesList += std::to_string(perGradeStartingIndex[i]);
        indicesList += (i==perGradeStartingIndex.size()-1)?"}":",";
    }

    return indicesList;
}


// build the string of the 'vectorNumber' basis name, prefixed by 'e' or 'E' according to 'lowercase'.
// example:  "E12", "e23", ...
std::string vectorIdToString(const MetaData& metaData, const bool lowercase, const int vectorNumber){
    return lowercase?"e"+metaData.basisVectorName[0]:"E"+metaData.basisVectorName[vectorNumber];
}


// generate a string that contains at table that summarize the result of the inner product between vectors (the metric /quadratic form of the algebra)
// example (c3ga):
// static const std::string metric =
// "    e0  e1  e2  e3  ei
// e0	0	0	0	0	-1
// e1	0	1	0	0	0
// e2	0	0	1	0	0
// e3	0	0	0	1	0
// ei	-1	0	0	0	0  ";
std::string metricToString(const MetaData &metaData){

    std::string metric = "\"\\\n\t";

    // first row : basis vectors
    for(auto &i: metaData.basisVectorName)
        metric += "e" + i + "\t";
    metric += "\\n\\\n"; // generates "\n\"

    // inner product table for vector basis
    for(unsigned int i=0; i<(unsigned int)metaData.metric.rows(); ++i){

        // start a table line by the vector name
        metric += "e" + metaData.basisVectorName[i] + "\t";

        for(unsigned int j=0; j<(unsigned int)metaData.metric.cols(); ++j){

            // if integer, remove the floating point notation
            if(metaData.metric(i,j) - (int)metaData.metric(i,j) == 0 )
                metric += std::to_string((int)metaData.metric(i,j)) + "\t";
            else metric += std::to_string(metaData.metric(i,j)) + "\t";
        }

        metric += "\\n\\\n"; // generates "\n\"
    }

    metric += "\"";

    return metric;
}


// generate an array containing the 'dimension+1' first binomial coefficients.
// This array is used in Mvec.hpp to know the size of a new k-vector to be added to the multivector.
// example (c3ga): static constexpr unsigned int binomialArray[6] = {1,5,10,10,5,1};
std::string binomialCoefToString(const unsigned int &dimension) {

    std::string tabBinomCoef="binomialArray["+std::to_string(dimension+1)+"] = "+"{";
    for(unsigned int i=0; i<=dimension; ++i){
        tabBinomCoef += std::to_string(bin_coeff(dimension,i));
        tabBinomCoef += (i==dimension)?"}":",";
    }

    return tabBinomCoef;
}



/*******************************************************************************************************
 ******************************** DUAL AND BASIS CHANGES HANDLING *************************************
 *******************************************************************************************************/


// generate a string containing the template of the function that fill the transformation matrices (direct here) for each grade of the vector space.
// This string is loaded in BasisTransform.hpp.
std::string initKVectorDirectTransformationMatrixFunctionPrototype(const unsigned int grade){
    // first generate the comments for the considered grade. 
    std::string outputString = "\t/// contains and load the components of the grade "+std::to_string(grade)+" direct transformation matrix \n ";
    // generate the prototype of the function that contains and load the components of the direct transformation matrix.
    outputString += "\ttemplate<typename T>\n";
    outputString +="\tEigen::SparseMatrix<T, Eigen::ColMajor>  loadgrade"+ std::to_string(grade)+"Matrix() {\n";
    outputString += "\t\tproject_grade"+std::to_string(grade)+"_basischange_components;\n";
    outputString += "\t\tstd::vector<T> gradeVectorComponents = decodeStringToVecOfT<T>(grade"+std::to_string(grade)+"MatrixComponents);\n";
    outputString += "\t\tproject_grade"+std::to_string(grade)+"_fill_matrix;\n"
                    "\t\treturn perGradeBasisTransformMatrix;\n"
                    "\t}\n";
    return outputString;
}


// generate a string containing the template of the function that fill the transformation matrices (inverse here) for each grade of the vector space.
// This string is loaded in BasisTransform.hpp.
std::string initKVectorInverseTransformationMatrixFunctionPrototype(const unsigned int grade){
    // first generate the comments for the considered grade. 
    std::string outputString = "\t/// contains and load the components of the grade "+std::to_string(grade)+" inverse transformation matrix \n ";
    // generate the prototype of the function that contains and load the components of the inverse transformation matrix.
    outputString += "\ttemplate<typename T>\n";
    outputString +="\tEigen::SparseMatrix<T, Eigen::ColMajor>  loadgrade"+ std::to_string(grade)+"InverseMatrix() {\n";
    outputString += "\t\tproject_grade"+std::to_string(grade)+"_basischange_components;\n";
    outputString += "\t\tstd::vector<T> gradeVectorComponents = decodeStringToVecOfT<T>(grade"+std::to_string(grade)+"MatrixComponents);\n";
    outputString += "\t\tproject_grade"+std::to_string(grade)+"_fill_matrix;\n"
                    "\t\treturn perGradeBasisTransformMatrix;\n"
                    "\t}\n";
    return outputString;
}




// generate a string containing the components of one transformation matrix (direct or inverse).
// This string is loaded in BasisTransform.hpp.
std::string perGradetransformMatricesToString(const std::vector<double>& transformComponents, const unsigned int grade,  const unsigned int startingPosition, const unsigned int endPosition) {

    std::string outputString = "const std::string grade"+std::to_string(grade)+"MatrixComponents = \" ";
    for(unsigned int j=startingPosition;j<endPosition;++j){
        outputString += std::to_string((int)transformComponents[3*j]) + " " + std::to_string((int)transformComponents[3*j+1]) + " " + std::to_string(transformComponents[3*j+2])+ " ";
    }
    outputString += "\""; // filling of the grade i components string is finished
    
    return outputString;
}

// generate a string containing the initialization of one transformation matrix (direct or inverse) using an eigen sparse vector.
// This string is loaded in BasisTransform.hpp.
std::string perGradeEigenSparseMatrixInitialization(const std::vector<double>& transformComponents, const unsigned int dimension, const unsigned int grade,  const unsigned int startingPosition, const unsigned int endPosition, const unsigned int numberOfNonNullComponents, const std::string& nameComponentsDataStructure) {

    std::string currentSize = std::to_string(bin_coeff(dimension,grade));
    std::string nameMatrixDataStructure =  "perGradeBasisTransformMatrix";
    std::string outputString="Eigen::SparseMatrix<T, Eigen::ColMajor> "+nameMatrixDataStructure;
    outputString += "\t\t = Eigen::SparseMatrix<T, Eigen::ColMajor>("+currentSize+","+currentSize+");\n";
    outputString += "\t\t"+ nameMatrixDataStructure + ".reserve(" + std::to_string(numberOfNonNullComponents) + ");\n";
    if(endPosition-startingPosition==1){
        outputString += "\t\t"+nameMatrixDataStructure + ".insert((int)"+nameComponentsDataStructure+"[0],(int)"+nameComponentsDataStructure+"[1]) = "+nameComponentsDataStructure+"[2]";
        return outputString;
    }
    outputString += "\t\tfor(unsigned int i=0;i<" + std::to_string(endPosition-startingPosition)+";++i)\n";
    outputString += "\t\t\t"+nameMatrixDataStructure + ".insert((int)"+nameComponentsDataStructure+"[3*i],(int)"+nameComponentsDataStructure+"[(3*i)+1]) = "+nameComponentsDataStructure+"[(3*i)+2]";
    return outputString;
}



// generate a string containing the initialization, filling of the transformation matrices either direct or inverse with respect to the boolean isInverse.
// This string is loaded in BasisTransform.hpp.
std::string loadAllDirectOrInverseMatrices(const std::vector<unsigned int> &sizeTransformationMatrices,
                                           const std::vector<double> &transformComponents, const bool isInverse){
    std::string outputString = "";
    unsigned int startingPosition = 0;
    unsigned int endPosition = 0; // enables to loop over the list of components

    // for each grade generate a string
    for(unsigned int i = 0; i < sizeTransformationMatrices.size();++i){
        endPosition += sizeTransformationMatrices[i];
        // generate the prototype of the loading of the current matrix
        std::string currentMatrixLoadingFunction = isInverse?initKVectorInverseTransformationMatrixFunctionPrototype(i):initKVectorDirectTransformationMatrixFunctionPrototype(i);
        // initialize the string containing the components
        substitute(currentMatrixLoadingFunction,"project_grade"+std::to_string(i)+"_basischange_components",perGradetransformMatricesToString(transformComponents, i,  startingPosition, endPosition));
        // now it is time to initialize the eigen matrix
        substitute(currentMatrixLoadingFunction,"project_grade"+std::to_string(i)+"_fill_matrix",perGradeEigenSparseMatrixInitialization(transformComponents, (unsigned int)(sizeTransformationMatrices.size()-1), i, startingPosition, endPosition, sizeTransformationMatrices[i], std::string("gradeVectorComponents")));
        startingPosition += sizeTransformationMatrices[i];

        outputString +=currentMatrixLoadingFunction+"\n";
    }
    outputString +="\n";

    return outputString;
}


// generate a string containing the initialization, filling of the transformation matrices either direct or inverse with respect to the boolean isInverse.
// This string is loaded in BasisTransform.hpp.
std::string callAllDirectOrInverseMatricesFunctions(const std::vector<unsigned int>& sizeTransformationMatrices, const bool isInverse){
    std::string outputString = "\ttemplate<typename T>\n";

    // the prototype of the function
    outputString +="\tconst std::array<Eigen::SparseMatrix<T, Eigen::ColMajor>,"+ std::to_string(sizeTransformationMatrices.size())+"> loadMatrices" + (isInverse?"Inverse()":"()") + " {\n";
    outputString +="\t\tstd::array<Eigen::SparseMatrix<T, Eigen::ColMajor>,"+ std::to_string(sizeTransformationMatrices.size())+">  transformationMatrices;\n";
                    

    // for each grade generate a string
    for(unsigned int grade = 0; grade < sizeTransformationMatrices.size();++grade){
        outputString += "\t\ttransformationMatrices["+std::to_string(grade)+"] = " +"load"+"grade"+std::to_string(grade)+(isInverse?"Inverse":"")+"Matrix<T>();\n";
    }

    outputString += "\t\treturn transformationMatrices;\n";
    outputString += "\t}\n";
    return outputString;
}


// generate a string containing the initialization, filling of the dual coefficients array 
// This string is loaded in DualCoefficients.hpp.
std::string loadAllDualCoefficientsArray(const std::vector<int> & sizeDualCoefficientsArray, const std::string& stringDualCoefficientsComponents){
    std::string outputString = "\ttemplate<typename T>\n";
    unsigned int startingPosition = 0;
    unsigned int endPosition = 0; // enables to loop over the list of components

    // the prototype of the function
    outputString +="\tconst std::array<Eigen::Matrix<T, Eigen::Dynamic,1>,"+ std::to_string(sizeDualCoefficientsArray.size())+"> loadFastDualArray() {\n";
    outputString +="\t\tstd::array<Eigen::Matrix<T, Eigen::Dynamic,1>,"+ std::to_string(sizeDualCoefficientsArray.size())+"> dualArrayCoefficients;\n";

    // put the string containing the coefficients of the dual
    outputString += "\t\tstd::string stringDualCoefficients=\""+stringDualCoefficientsComponents+"\";\n";
    outputString += "\t\tconst std::vector<T> vectorDualComponents = decodeStringOfDualCoefToVecOfT<T>(stringDualCoefficients);\n";

    std::string nameComponentsDataStructure = "vectorDualComponents";
    std::string nameArrayDataStructure = "dualArrayCoefficients";

    // for each grade generate a string
    for(unsigned int i = 0; i < sizeDualCoefficientsArray.size()-1;++i){
        startingPosition = sizeDualCoefficientsArray[i];
        endPosition = sizeDualCoefficientsArray[i+1];
        std::string currentSize = std::to_string(endPosition-startingPosition);
        outputString += "\t\t" +nameArrayDataStructure  + "[" + std::to_string(i) + "]= Eigen::Matrix<T, Eigen::Dynamic,1>("+currentSize+");\n";

        outputString += "\t\tfor(unsigned int i=" + std::to_string(startingPosition) + ";i<" + std::to_string(endPosition)+";++i){\n";
        outputString += "\t\t\t"+nameArrayDataStructure  + "[" + std::to_string(i) + "].coeffRef(i-"+ std::to_string(startingPosition) + ") = (T)"+nameComponentsDataStructure+"[i]; \n";
        outputString += "\t\t}\n";
    }
    // Last grade
    startingPosition = sizeDualCoefficientsArray.back();
    endPosition = startingPosition+1;
    std::string currentSize = std::to_string(endPosition-startingPosition);
    outputString += "\t\t" +nameArrayDataStructure  + "[" + std::to_string(sizeDualCoefficientsArray.size()-1) + "]= Eigen::Matrix<T, Eigen::Dynamic,1>("+currentSize+");\n";

    outputString += "\t\tfor(unsigned int i=" + std::to_string(startingPosition) + ";i<" + std::to_string(endPosition)+";++i){\n";
    outputString += "\t\t\t"+nameArrayDataStructure  + "[" + std::to_string(sizeDualCoefficientsArray.size()-1) + "].coeffRef(0) = (T)"+nameComponentsDataStructure+"[i]; \n";
    outputString += "\t\t}\n";
    outputString += "\t\treturn dualArrayCoefficients;\n";
    outputString +="\t}\n";

    return outputString;
}




//////////////////////////////////////////////////////// END FAST DUAL COEFFICIENTS AND TRANSFORMATION MATRICES //////////////////////////////////////////////////////////////



// prototype of the methods namespace::e12() in Python
std::string staticOneComponentMultivectorPrototypePython(){
    std::string res = "m.def(\"eproject_name_blade\", &eproject_name_blade<double>);\n";
    return res;
}




// prototype of the methods namespace::e12()
// the purpose is to generate the methods designed to assign to a specific blade of the multivector a component (scalar).
// example (c3ga): Mvec a = 3 * namespace::e12();
std::string staticOneComponentMultivectorPrototype(){
    std::string res = "    /// \\brief return a multivector that contains only the unit basis k-vector project_name_blade.\n";
    res += "    template<typename T>\n    static Mvec<T> eproject_name_blade(){\n        Mvec<T> res;\n        return res.componentToOne(project_grade_blade, project_homogeneous_index_blade);\n    }\n\n";
    return res;
}

// prototype of the methods mv.e12()
std::string oneComponentMultivectorPrototype(){
    std::string res = "        Mvec eproject_name_blade() const {return this->extractOneComponent(project_grade_blade,project_current_binomial_coefficient, project_homogeneous_index_blade);}\n";
    return res;
}

// definitions of the constants
// the purpose is to generate constants values defining the index of each basis vector in the multivector (from 0 for scalar to 2^d-1 for pseudoscalar)
// example (c3ga) :
// const unsigned int scalar = 0;
// const unsigned int E0  = 1;
// const unsigned int E1  = 2;
// const unsigned int E2  = 4;
// const unsigned int E12 = 6;
// ...
std::string constantsDefinition(){
    std::string res = "    const unsigned int Eproject_name_blade = project_xor_index_blade;\n";
    return res;
}

// in case the metric is not full rank, the pseudo scalar will be 0, thus it is not possible to compute some operation like dual.
// this functions will thus be commented on the source code
std::string singularMetricCommentBegin(){
    std::string res = "    /* the metric defined in this algebra is not full ranked, so the following fuction can not be implemented";
    return res;
}

// in case the metric is not full rank, the pseudo scalar will be 0, thus it is not possible to compute some operation like dual.
// this functions will thus be commented on the source code
std::string singularMetricCommentEnd(){
    std::string res = "    */\n";
    return res;
}

// in case the metric is not full rank, the pseudo scalar will be 0, thus it is not possible to compute some operation like dual.
// this functions will thus be commented on the source code

// compute the sequence of signs used in the reverse operation
// if mv = a1^a2^...^an, then reverse(mv) = an^...^a2^a1 = (-1)^(0.5*grade*(grade-1)) mv
// These grade dependent signs are stored in an array
// example (c3ga): static const int signReversePerGrade[6]={1,1,-1,-1,1,1};
std::string reverseSignArrayToString(const unsigned int &dimension){

    std::string resSignReverse = "signReversePerGrade["+std::to_string(dimension+1)+"] = {";

    for(unsigned int i=0;i<=dimension;++i){
        resSignReverse += std::to_string((int)pow(-1,(i*(i-1))/2));
        if(i<dimension) resSignReverse += ",";
    }

    return resSignReverse + "}";
}


// convert to string the array associated the xor indices to the grade and the index in the VectorXd of the multivector
std::string xorIndexToGradeAndHomogeneousIndexArraysToString(unsigned int dimension, const ProductTools& product){
    std::string outputString = "";

    std::string stringXorToGrade="    constexpr unsigned int xorIndexToGrade[] = {";
    std::string stringXorToHomogeneousIndex="    constexpr unsigned int xorIndexToHomogeneousIndex[] = {";

    // for each xor index in the range 0 to 2^dim
    for(unsigned int xorIndices=0;xorIndices<(unsigned int)((1<<dimension)-1);++xorIndices){
        // add the grade and the homogeneous index to their corresponding array
        stringXorToGrade += std::to_string(product.getGrade(xorIndices)) + ",";
        stringXorToHomogeneousIndex += std::to_string(product.getHomogeneousIndex(xorIndices)) + ",";
    }
    // close the two arrays and comment the two data structures
    stringXorToGrade+= std::to_string(product.getGrade((unsigned int)((1<<dimension)-1))) + "}; /*!< given a Xor index in a multivector, this array indicates the corresponding grade*/ \n\n";
    stringXorToHomogeneousIndex+= std::to_string(product.getHomogeneousIndex((unsigned int)((1<<dimension)-1))) + "}; /*!< given a Xor index in a multivector, this array indicates the corresponding index in the whole homogeneous vector*/";


    outputString += stringXorToGrade + stringXorToHomogeneousIndex;
    return outputString;
}

// defines the array of permutations and coefficients required in the computation of the fast dual
std::string fastDualUtilitiesBasisChange(unsigned int dimension, const ProductTools& product,
                                         const std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor> >& transformationMatrices,
                                         const std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor> >& inverseTransformationMatrices,
                                         const Eigen::VectorXd &diagonalMetric,
                                         double scaleInversePseudoScalar,
                                         std::string srcDirectory,
                                         std::string& fastDualComponents){
    std::string outputStringPermutations="    const std::array<std::vector<unsigned int>, " + std::to_string(dimension+1) + "> dualPermutations = {{ ";

    std::string outputStringCoefficients="\n    const std::array<Eigen::Matrix<double, Eigen::Dynamic,1>, "+std::to_string(dimension+1)+"> dualCoefficients = loadFastDualArray<double>(); /*!< array containing some basis change coefficients required to compute the dual */\n    "; // replace by load

    std::vector<double> dualCoefficientsComponents;


    //std::string outputStringCoefficients="\n\n    template<typename T>\n"
    //        "    Eigen::Matrix<T, Eigen::Dynamic, 1> coefficientsVector = []()->Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>{Eigen::Matrix<T, Eigen::Dynamic, 1> tmp("+std::to_string(1<<dimension)+");tmp<<";
    // for each possible grade from 0 to dimension-1 do
    for(unsigned int grade=0;grade<dimension;++grade){
        if(grade==0) outputStringPermutations += "{ ";
        else outputStringPermutations += "{{ ";

        //outputStringCoefficients+="[]()->Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>{Eigen::Matrix<T, Eigen::Dynamic, 1> tmp("+std::to_string(bin_coeff(dimension,grade))+");tmp<<";

        std::list<productComponent<double>> listProductInnerProduct =product.generateExplicitInnerProductList(
                grade,
                dimension,
                inverseTransformationMatrices[grade],
                inverseTransformationMatrices[dimension],
                transformationMatrices[dimension-grade],
                diagonalMetric);


        listProductInnerProduct.sort([](const productComponent<double> &pc1, const productComponent<double> &pc2) {return pc2.indexOfMv1 > pc1.indexOfMv1;});
        for(auto & products : listProductInnerProduct){ // for each list of permutation for the grade 'grade'
            outputStringPermutations += std::to_string(products.indexOfMv3);
            if(&products != &listProductInnerProduct.back()){
                outputStringPermutations += ",";// add the permutation required to compute the dual
            }else{
                if(grade==0) outputStringPermutations += "}, ";
                else outputStringPermutations += "}}, ";

            }
        }

        listProductInnerProduct.sort([](const productComponent<double> &pc1, const productComponent<double> &pc2) {return pc2.indexOfMv3 > pc1.indexOfMv3;});
        for(auto & products : listProductInnerProduct){ // for each list of permutation for the grade 'grade'
            //outputStringCoefficients += std::to_string(products.coefficient);
            dualCoefficientsComponents.push_back(products.coefficient);
        }
    }

    // last element correspond to permutation required for the pseudo scalar (grade = dimension)
    std::list<productComponent<double>> listProductInnerProduct =product.generateExplicitInnerProductList(
            dimension,
            dimension,
            inverseTransformationMatrices[dimension],
            inverseTransformationMatrices[dimension],
            transformationMatrices[0],
            diagonalMetric);

    //outputStringCoefficients += "[]()->Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>{Eigen::Matrix<T, Eigen::Dynamic, 1> tmp(1);tmp<<";
    outputStringPermutations += "{";
    listProductInnerProduct.sort([](const productComponent<double> &pc1, const productComponent<double> &pc2) {return pc2.indexOfMv1 > pc1.indexOfMv1;});


    // when the metric is singular, there might be no components at all then this loop might be ignored
    // to avoid this, we put one element with a 0 into the two structures.
    if(listProductInnerProduct.empty()){
        outputStringPermutations += std::to_string(0)+"} }};\n";
        dualCoefficientsComponents.push_back(0);
    }  /// \todo should not happen anymore

    for(auto& products : listProductInnerProduct){ // for each list of permutation for the grade 'grade'
        outputStringPermutations += std::to_string(products.indexOfMv3)+"} }}; /*!< array referring to some permutations required to compute the dual. */\n";// add the permutation required to compute the dual
        dualCoefficientsComponents.push_back(products.coefficient);
    }

    for(unsigned int i=0;i<dualCoefficientsComponents.size();++i){
        fastDualComponents+= std::to_string(dualCoefficientsComponents[i])+" ";
    }


    return outputStringPermutations + outputStringCoefficients;
}



// defines the array of permutations and coefficients required in the computation of the fast dual
// reproduce the same behavior than with basis change
std::string fastDualUtilities(unsigned int dimension, const ProductTools& product,
                              const Eigen::VectorXd &diagonalMetric,
                              double scaleInversePseudoScalar,
                              std::string srcDirectory,
                              std::string& fastDualComponents){

    std::string outputStringPermutations="    std::array<std::vector<unsigned int>, " + std::to_string(dimension+1) + "> dualPermutations = {{ ";

    // this string was previously initialized at library compile time. Now this is done in a similar way to the construction of the transformation matrices
    // loading it in a file
    std::string outputStringCoefficients = "\n"
                                           "    const std::array<Eigen::Matrix<double, Eigen::Dynamic,1>, " + std::to_string(dimension+1) + "> dualCoefficients = loadFastDualArray<double>(); /*!< array containing the coefficients needed to compute the dual */ \n    ";

    std::vector<double> dualCoefficientsComponents;


    //std::string outputStringCoefficients="\n\n    template<typename T>\n"
    //        "    Eigen::Matrix<T, Eigen::Dynamic, 1> coefficientsVector = []()->Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>{Eigen::Matrix<T, Eigen::Dynamic, 1> tmp("+std::to_string(1<<dimension)+");tmp<<";
    // for each possible grade from 0 to dimension-1 do
    for(unsigned int grade=0;grade<dimension;++grade){
        if(grade==0) outputStringPermutations += "{ ";
        else outputStringPermutations += "{{ ";

        //outputStringCoefficients+="[]()->Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>{Eigen::Matrix<T, Eigen::Dynamic, 1> tmp("+std::to_string(bin_coeff(dimension,grade))+");tmp<<";
        std::list<productComponent<double>> listProductInnerProduct =product.generateExplicitInnerProductListEuclideanSpace(
                grade,
                dimension,
                diagonalMetric);
        
        // Compute the permutation required 
        listProductInnerProduct.sort([](const productComponent<double> &pc1, const productComponent<double> &pc2) {return pc2.indexOfMv1 > pc1.indexOfMv1;});
        for(auto & products : listProductInnerProduct){ // for each list of permutation for the grade 'grade'
            outputStringPermutations += std::to_string(products.indexOfMv3);
            if(&products != &listProductInnerProduct.back()){
                outputStringPermutations += ", ";// add the permutation required to compute the dual
            }else{
                if(grade==0) outputStringPermutations += "}, ";
                else outputStringPermutations += "}}, ";
            }
        }

        // Compute the coefficient required 
        listProductInnerProduct.sort([](const productComponent<double> &pc1, const productComponent<double> &pc2) {return pc2.indexOfMv3 > pc1.indexOfMv3;});
        for(auto & products : listProductInnerProduct){ // for each list of permutation for the grade 'grade'
            //outputStringCoefficients += std::to_string(products.coefficient);
            dualCoefficientsComponents.push_back(products.coefficient);

        }

    }
    // last element correspond to permutation required for the pseudo scalar (grade = dimension)
    std::list<productComponent<double>> listProductInnerProduct =product.generateExplicitInnerProductListEuclideanSpace(
            dimension,
            dimension,
            diagonalMetric);

    //outputStringCoefficients += "[]()->Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>{Eigen::Matrix<T, Eigen::Dynamic, 1> tmp(1);tmp<<";
    outputStringPermutations += "{";
    listProductInnerProduct.sort([](const productComponent<double> &pc1, const productComponent<double> &pc2) {return pc2.indexOfMv1 > pc1.indexOfMv1;});

    // when the metric is singular, there might be no components at all then this loop might be ignored
    // to avoid this, we put one element with a 0 into the two structures.
    if(listProductInnerProduct.empty()){
        outputStringPermutations += std::to_string(0)+"} }};\n";
        dualCoefficientsComponents.push_back(0);
    } /// \todo should not hapen anymore

    for(auto& products : listProductInnerProduct){ // for each list of permutation for the grade 'grade'
        outputStringPermutations += std::to_string(products.indexOfMv3)+"} }};\n";// add the permutation required to compute the dual
        //outputStringCoefficients += std::to_string(products.coefficient)+";return tmp;}()\n    ";// add the coefficient required to compute the dual
        dualCoefficientsComponents.push_back(products.coefficient);
    }
    //outputStringCoefficients += "}};";
    
    for(unsigned int i=0;i<dualCoefficientsComponents.size();++i){
        fastDualComponents+= std::to_string(dualCoefficientsComponents[i])+" ";
    }



    return outputStringPermutations + outputStringCoefficients;
}







// defines the array of permutations and coefficients required in the computation of the fast dual
std::string primalWedgeDualUtilities(const unsigned int dimension, const ProductTools& product,
                                     const Eigen::VectorXd &diagonalMetric,
                                     double scaleInversePseudoScalar){

    std::string outputString = "";
    outputString += "    template<typename T>\n"
                    "    std::array<T, " + std::to_string(1<<dimension) + "> recursiveDualCoefficients = {{ ";

    std::vector<double> tabDualCoefficients(1<<dimension);

    // for each possible grade from 0 to dimension-1 do
    for(unsigned int grade=0; grade<=dimension; ++grade){

        std::list<productComponent<double>> listProductInnerProduct = product.generateExplicitInnerProductListEuclideanSpace(
                grade,
                dimension,
                diagonalMetric);

        for(auto & products : listProductInnerProduct) // for each list of permutation for the grade 'grade'
            tabDualCoefficients[product.getXorIndex(grade,products.indexOfMv3)] = products.coefficient;
    }

    for(int i=0 ; i<(1<<dimension)-1; ++i){ // for each list of permutation for the grade 'grade'
        outputString += std::to_string(tabDualCoefficients[i]);
        outputString += ", ";// add the coefficient required to compute the dual
    }
    outputString += std::to_string(tabDualCoefficients[0]) +"}}; /*!< array containing the coefficients needed to compute the recursive product like (primal^dual) */\n    ";
    return outputString;
}


// defines the array of permutations and coefficients required in the computation of the fast dual
std::string primalWedgeDualUtilitiesBasisChange(unsigned int dimension, const ProductTools& product,
                                                const std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor> >& transformationMatrices,
                                                const std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor> >& inverseTransformationMatrices,
                                                const Eigen::VectorXd &diagonalMetric,
                                                double scaleInversePseudoScalar){
    std::string outputString = "    template<typename T>\n"
                               "    std::array<T, " + std::to_string(1<<dimension) + "> recursiveDualCoefficients = {{ ";

    std::vector<double> tabDualCoefficients(1<<dimension);

    // for each possible grade from 0 to dimension-1 do
    for(unsigned int grade=0; grade<=dimension; ++grade){

        std::list<productComponent<double>> listProductInnerProduct = product.generateExplicitInnerProductList(
                grade,
                dimension,
                inverseTransformationMatrices[grade],
                inverseTransformationMatrices[dimension],
                transformationMatrices[dimension-grade],
                diagonalMetric);

        for(auto & products : listProductInnerProduct){ // for each list of permutation for the grade 'grade'
            tabDualCoefficients[product.getXorIndex(grade,products.indexOfMv3)] = products.coefficient;
        }
    }

    for(int i=0 ; i<(1<<dimension)-1; ++i){ // for each list of permutation for the grade 'grade'
        outputString += std::to_string(tabDualCoefficients[i]);
        outputString += ",";// add the coefficient required to compute the dual
    }
    outputString += std::to_string(tabDualCoefficients[0]) +"}}; /*!< array containing the coefficients needed to compute the recursive product like (primal^dual) */";

    return outputString;
}








// generate the prototype of BasisTransformComponents
std::string basisTransformMatricesLoad(){
    std::string outputString="";
    outputString += "    template<typename T>\n    const std::array<Eigen::SparseMatrix<T, Eigen::ColMajor>,project_dim_plus_one> transformationMatrices = loadMatrices<T>(); /*!< set transformation matrices to transform a k-vector from the orhogonal basis to the original basis */\n";
    outputString += "    template<typename T>\n    const std::array<Eigen::SparseMatrix<T, Eigen::ColMajor>,project_dim_plus_one> transformationMatricesInverse = loadMatricesInverse<T>(); /*!< set transformation matrices to transform a k-vector from the original basis to the orhogonal basis */\n";

    return outputString;
}


// generate the methods and constants designed to build multivectors that contain only one specified component of the current multivector.
// More precisely, all constant blades identifier (e.g 12, 13, 14 ...),  their corresponding position in both the general multivector (XOR index) and the homogeneous multivector
// this function enables to initialize the static and non static functions e12() for example,
// as well as the constants blades, for example const unsigned E12 = 3
std::string multivectorComponentBuilder(const MetaData& metaData, const std::string &data){

    // output generated data
    std::string outputString = {};

    // for each grade from 1 to the dimension (the scalar is already generated)
    for(unsigned int grade=1; grade<=metaData.dimension; ++grade) {

        // for high dimension, do not generate all e(...)
        if(bin_coeff(metaData.dimension,grade) > metaData.maxDimBasisAccessor)
            continue;

        // for the current grade, generate all the XOR indices
        std::vector<bool> booleanXorIndex(metaData.dimension);
        std::fill(booleanXorIndex.begin(), booleanXorIndex.end(), false);

        // build the first xorIndex of grade 'grade' (with 'grade' values to true).
        std::fill(booleanXorIndex.begin(), booleanXorIndex.begin() + grade, true);

        // for each permutation of the true values of booleanXorIndex
        unsigned int positionInKVector = 0;
        do {
            // initialize a template string
            std::string currentData = data;

            // convert the vector of bool to a unsigned int (each bit get a boolean value) => xor index
            unsigned int xorIndex = 0;
            std::string basisBlades = ""; // example : e13
            for(unsigned int i=0; i<metaData.dimension; ++i) {
                if(booleanXorIndex[i]) {
                    xorIndex += (1<<i); // Xor index
                    basisBlades += metaData.basisVectorName[i];
                }
            }

            // edit the template string
            substitute(currentData,"project_name_blade", basisBlades);
            substitute(currentData,"project_homogeneous_index_blade", std::to_string(positionInKVector));
            substitute(currentData,"project_xor_index_blade", std::to_string(xorIndex));
            substitute(currentData,"project_grade_blade", std::to_string(grade));
            substitute(currentData,"project_current_binomial_coefficient", std::to_string(bin_coeff(metaData.dimension,grade)));
            //basisBladesStr += " = " + std::to_string(xorIndex) + ";\n";
            positionInKVector++;

            // push back the last block of code where we substitute some parameters
            outputString += currentData;

        } while(std::prev_permutation(booleanXorIndex.begin(), booleanXorIndex.end())); // compute next permutation of the true values of booleanXorIndex
    }

    return outputString;
}




// the diagonal metric is stored as a vector
std::string diagonalMetricToString(const MetaData& metaData){

    // vector definition
    std::string outputString = "template<typename T>\n    Eigen::Matrix<T, Eigen::Dynamic, 1> diagonalMetric = []()->Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>{Eigen::Matrix<T, Eigen::Dynamic, 1> tmp("+std::to_string(metaData.dimension)+");tmp<<";

    // add vector elements
    for(unsigned int i=0; i<(unsigned int)metaData.diagonalMetric.size()-1; ++i)
        outputString += std::to_string(metaData.diagonalMetric(i)) + ",";

    outputString += std::to_string(metaData.diagonalMetric(metaData.diagonalMetric.size()-1)) + "; return tmp;}();";

    // add doxygen comments
    outputString += "   /*!< defines the diagonal metric (stored as a vector) */";

    return outputString;
}


/// \todo TO REMOVE change the generated function from "project" to "toOrthogonalBasis" and "toOriginalBasis"
std::string templateVectorToOrthogonalBasisToString(){
    std::string result = "Mvec project() const {\n"
            "            Mvec mvProjected;\n"
            "            for(auto & itMv : this->mvData)\n"
            "                mvProjected.mvData[itMv.first] = (transformationMatricesInverse<T>[itMv.first]*itMv.second); // (transformationMatricesInverse<T>[itMv.first]*itMv.second)\n"
            "            return mvProjected;\n"
            "        }\n"
            "\n"
            "        Mvec backProject() const {\n"
            "            Mvec mvBackProjected;\n"
            "            for(auto & itMv : this->mvData){\n"
            "                if((itMv.second.array() != 0.0).any()){\n"
            "                    // Compute the change of basis for mv2: from non-orthogonal to orthogonal\n"
            "                    mvBackProjected.mvData[itMv.first] = (transformationMatrices<T>[itMv.first] * itMv.second);\n"
            "                }\n"
            "            }\n"
            "            return mvBackProjected;\n"
            "        }";
    return result;
}

std::string templateTransformVectorToOrthogonalBasisToString(){
    std::string result = "Mvec<T> mv1Ortho;\n"
            "        Mvec<T> mv2Ortho;\n"
            "        Mvec<T> resultOrtho;\n"
            "\n"
            "        // transformation matrix applied to mv1\n"
            "        mv1Ortho = this->project();\n"
            "        mv2Ortho = mv2.project();";
    return result;
}

std::string templateTransformVectorToOriginalBasisToString(){
    std::string result = "mv3 = resultOrtho.backProject();";
    return result;
}




/*******************************************************************************************************
 ******************************** EXPLICIT OUTER PRODUCT DEFINITION ************************************
 *******************************************************************************************************/

std::string outerProductExplicitComments(const unsigned int gradeMv1, const unsigned int gradeMv2, const unsigned int gradeMv3) {
    std::string outerFunctionComments = "/// \\brief Compute the outer product between two homogeneous multivectors mv1 (grade "+ std::to_string(gradeMv1)+") and mv2 (grade " + std::to_string(gradeMv2) + "). \n\t";
    outerFunctionComments += "/// \\tparam the type of value that we manipulate, either float or double or something.\n\t";
    outerFunctionComments += "/// \\param mv1 - the first homogeneous multivector of grade " + std::to_string(gradeMv1) + " represented as an Eigen::VectorXd\n\t";
    outerFunctionComments += "/// \\param mv2 - the second homogeneous multivector of grade " + std::to_string(gradeMv2) + " represented as a Eigen::VectorXd\n\t";
    outerFunctionComments += "/// \\param mv3 - the result of mv1^mv2, which is also a homogeneous multivector of grade " + std::to_string(gradeMv3) + "\n\t";
    return outerFunctionComments;
}


std::string outerProductExplicitPrototype(const unsigned int gradeMv1, const unsigned int gradeMv2){
    return "template<typename T>\n\tvoid outer_"+std::to_string(gradeMv1)+"_"+std::to_string(gradeMv2)+"(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){\n";
}


std::string productListToString(std::list<productComponent<double>> &listOfExplicitProducts) {

    // declare an empty string (to use += later on)
    std::string outputString = {};

    // for a C++ generated function, the list is sorted according to the index of the result.
    listOfExplicitProducts.sort([](const productComponent<double> &pc1, const productComponent<double> &pc2) {return pc2.indexOfMv3 > pc1.indexOfMv3;});

    // for each component of mv3 = mv1 prod mv2
    auto iteratorProduct = listOfExplicitProducts.begin();
    while(iteratorProduct != listOfExplicitProducts.end()){

        // start an new component of mv3
        outputString += "\t\tmv3.coeffRef(" + std::to_string(((*iteratorProduct).indexOfMv3)) + ") +=";
        bool firstProduct = true;

        // compute how many product from mv1 and mv2 contribute to this current component of mv3
        productComponent<double> currentComponent = *iteratorProduct;
        unsigned long numberOfOccurrences = std::count_if(iteratorProduct, listOfExplicitProducts.end(),[currentComponent](productComponent<double> product){return product.indexOfMv3 == currentComponent.indexOfMv3;});

        // add all the products that contribute to this current component of mv3
        for(unsigned long i=0; i<numberOfOccurrences ; ++i){
            // put the coefficient
            // if the coefficient is either 1 or -1 then it can be associated to a sign
            if(firstProduct){
                // If the coefficient is neither -1 nor 1 then it has to appear in the generated component
                if( ((*iteratorProduct).coefficient !=-1) && ((*iteratorProduct).coefficient !=1)){
                    outputString += " "+std::to_string((*iteratorProduct).coefficient) + "*";
                }else{
                    outputString += ((*iteratorProduct).coefficient < 0) ? " -" : "  ";
                }
                firstProduct = false;
            }else {
                if (((*iteratorProduct).coefficient != -1) && ((*iteratorProduct).coefficient != 1)) {
                    outputString += ((*iteratorProduct).coefficient < 0) ? " " : " + ";
                    outputString += std::to_string((*iteratorProduct).coefficient) + "*";
                } else outputString += ((*iteratorProduct).coefficient < 0) ? " - " : " + ";
            }
            // write the product
            outputString += "mv1.coeff(" + std::to_string((*iteratorProduct).indexOfMv1) + ")" + "*mv2.coeff(" + std::to_string((*iteratorProduct).indexOfMv2)+ ")";

            // next product
            iteratorProduct++;
        }

        // next component of mv3, add a separator
        outputString += ";\n";
    }
    return outputString;
}


std::string generateOuterRecursive(const unsigned int gradeMv1, const unsigned int gradeMv2, const unsigned int gradeMv3){

    // call the fine recursive product
    return std::string("\t\touterProductHomogeneous(mv1,mv2,mv3," + std::to_string(gradeMv1) + "," + std::to_string(gradeMv2) + "," + std::to_string(gradeMv3) + ");\n");
}


// generate a string that contains all the c++ functions of the outer product
std::string generateOuterExplicit_cpp(const MetaData &metaData, const ProductTools& product){

    // string that will contain the result
    std::string outputString = "";

    // for all possible pair of grade (gradeMv1,gradeMv2)
    for(unsigned int gradeMv1=0; gradeMv1<=metaData.dimension; ++gradeMv1){
        for(unsigned int gradeMv2=0; gradeMv2<=metaData.dimension; ++gradeMv2){

            // discard out of grade outer product
            if((gradeMv1+gradeMv2) > metaData.dimension)
                continue;

            // current function: generate comments
            outputString += outerProductExplicitComments(gradeMv1, gradeMv2, gradeMv1+gradeMv2);

            // current function: generate prototype
            outputString += "template<typename T>\n\tvoid outer_"+ std::to_string(gradeMv1) + "_" + std::to_string(gradeMv2) + "(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){\n";

            // if Mv1 or Mv2 is a scalar, then use the eigen vector to scalar product
            if( (gradeMv1==0) ^ (gradeMv2==0) ){
                if(gradeMv1==0) outputString += "\t\tmv3 += mv1.coeff(0)*mv2;\n";
                else outputString += "\t\tmv3 += mv1*mv2.coeff(0);\n";
            } else { // not a scalar to mv product. Do a regular optimization

                // if one of the k-vector mv1 or mv2 is higher than the threshold, we use the recursive function.
                if( (bin_coeff(metaData.dimension,gradeMv1) > metaData.maxDimPrecomputedProducts) || (bin_coeff(metaData.dimension,gradeMv2) > metaData.maxDimPrecomputedProducts) ){

                    // generate a call to the fine recursive outer product
                    outputString += generateOuterRecursive(gradeMv1,gradeMv2,gradeMv1+gradeMv2);
                }else{

                    // compute all the components of the current product mv1{gradeMv1} ^ mv2{gradeMv2}.
                    // We get a list of quadruplet (indexMv1,indexMv2,indexMv3,coefficient),
                    // meaning Mv3{gradeMv3}[indexMv3] += coefficient * Mv1{gradeMv1}[indexMv1] * Mv2{gradeMv2}[indexMv2]
                    std::list<productComponent<double>> listOfProducts = product.generateExplicitOuterProductList(gradeMv1, gradeMv2);

                    // convert the products to c++ instructions
                    outputString += productListToString(listOfProducts);
                }
            }

            // end of the current function, finishing with a "}"
            outputString += "\t}\n\n\n\t";
        }
    }

    return outputString;
}


// generate an array of function pointers referring to each explicit products.
// Start with the declaration of the 2D array of functions pointer, then fill it.
// Note that the size of this functions pointer is (dimension+1)*(dimension+1) because a mv can have (dimension+1) kvectors
std::string generateOuterExplicitFunctionsPointer(const unsigned int &dimension){

    // string that will contain the result
    std::string outputString = "template<typename T>\n"
                                       "\tstd::array<std::array<std::function<void(const Eigen::Matrix<T, Eigen::Dynamic, 1> & , "
                                       "const Eigen::Matrix<T, Eigen::Dynamic, 1> & , Eigen::Matrix<T, Eigen::Dynamic, 1>&)>, "+ std::to_string(dimension+1)+">, "+ std::to_string(dimension+1)+"> outerFunctionsContainer = {{\n";

    // for all couple of grades
    for(unsigned int gradeMv1=0; gradeMv1<=dimension; ++gradeMv1){

        // beginning of a line
        outputString += "\t\t{{";

        for(unsigned int gradeMv2=0; gradeMv2<=dimension; ++gradeMv2){

            // when the grade of mv1^mv2 is higher than the dimension, we create a null pointer function
            if((gradeMv2+gradeMv1)>dimension) outputString += "{}";
                // else we create a pointer to the fine explicit outer product function which is identified using gradeMv1 and gradeMv2
            else outputString += "outer_" + std::to_string(gradeMv1) + "_" + std::to_string(gradeMv2) + "<T>";

            // Close the array when gradeMv2 reaches 'dimension'
            // else if gradeMv1 is 'dimension', then put a carriage return, else put a comma
            outputString += (gradeMv2<dimension)?",":"}}"+((gradeMv1==dimension)?std::string("\n"):std::string(",\n"));
        }
    }

    // when all the pointers are initialized, close the array of array with a double }
    outputString += "\t}};";

    return outputString;
}

















/*******************************************************************************************************
 ******************************** EXPLICIT INNER PRODUCT DEFINITION ************************************
 *******************************************************************************************************/

std::string innerProductExplicitComments(const unsigned int gradeMv1, const unsigned int gradeMv2, const unsigned int gradeMv3) {
    std::string outerFunctionComments = "/// \\brief Compute the inner product between two homogeneous multivectors mv1 (grade "+ std::to_string(gradeMv1)+") and mv2 (grade " + std::to_string(gradeMv2) + "). \n\t";
    outerFunctionComments += "/// \\tparam the type of value that we manipulate, either float or double or something.\n\t";
    outerFunctionComments += "/// \\param mv1 - the first homogeneous multivector of grade " + std::to_string(gradeMv1) + " represented as an Eigen::VectorXd\n\t";
    outerFunctionComments += "/// \\param mv2 - the second homogeneous multivector of grade " + std::to_string(gradeMv2) + " represented as a Eigen::VectorXd\n\t";
    outerFunctionComments += "/// \\param mv3 - the result of mv1.mv2, which is also a homogeneous multivector of grade " + std::to_string(gradeMv3) + "\n\t";
    return outerFunctionComments;
}

std::string generateInnerRecursiveBasisChange(const unsigned int gradeMv1, const unsigned int gradeMv2, const unsigned int gradeMv3,
                                              const unsigned int dimension) {

    //std::string outputString = "\t\tmv3 = Mvec<T>::transformationMatricesInverse[" + std::to_string(gradeMv3) + "]*mv3;\n";
    std::string outputString = "\t\tEigen::Matrix<T, Eigen::Dynamic, 1> tmp = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(" + std::to_string(bin_coeff(dimension,gradeMv3)) + ");\n";
    // call the fine recursive product (right or left contraction)
    if(gradeMv3 == (gradeMv2-gradeMv1))
        outputString += "\t\tleftContractionProductHomogeneous<T>";
    else
        outputString += "\t\trightContractionProductHomogeneous<T>";

    // put the arguments
    outputString += "(transformationMatricesInverse<T>[" + std::to_string(gradeMv1) + "]*mv1,"
                    + "transformationMatricesInverse<T>[" + std::to_string(gradeMv2) + "]*mv2,tmp,"
                    + std::to_string(gradeMv1) + "," + std::to_string(gradeMv2) + "," + std::to_string(gradeMv3) + ");\n";

    // change basis for mv3
    outputString += "\t\tmv3 += transformationMatrices<T>[" + std::to_string(gradeMv3) + "]*tmp;\n";

    return outputString;
}


std::string generateInnerRecursiveBasisChangeFloat64(const unsigned int gradeMv1, const unsigned int gradeMv2, const unsigned int gradeMv3, const unsigned int dimension) {

    //std::string outputString = "\t\tmv3 = Mvec<double>::transformationMatricesInverse[" + std::to_string(gradeMv3) + "]*mv3;\n";
    std::string outputString = "\t\tEigen::Matrix<double, Eigen::Dynamic, 1> tmp = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(" + std::to_string(bin_coeff(dimension,gradeMv3)) + ");\n";
    // call the fine recursive product (right or left contraction)

    if(gradeMv3 == (gradeMv2-gradeMv1))
        outputString += "\t\tleftContractionProductHomogeneous<double>";
    else
        outputString += "\t\trightContractionProductHomogeneous<double>";

    // put the arguments
    outputString += "(transformationMatricesInverse<T>[" + std::to_string(gradeMv1) + "]*mv1,"
                    + "transformationMatricesInverse<T>[" + std::to_string(gradeMv2) + "]*mv2,tmp,"
                    + std::to_string(gradeMv1) + "," + std::to_string(gradeMv2) + "," + std::to_string(gradeMv3) + ");\n";

    // change basis for mv3
    outputString += "\t\tmv3 += transformationMatrices<T>[" + std::to_string(gradeMv3) + "]*tmp;\n";

    return outputString;
}


std::string generateInnerRecursive(const unsigned int gradeMv1, const unsigned int gradeMv2, const unsigned int gradeMv3) {

    // call the fine recursive product (right or left contraction)
    if(gradeMv3 == (gradeMv2-gradeMv1))
        return std::string("\t\tleftContractionProductHomogeneous<T>(mv1,mv2,mv3," + std::to_string(gradeMv1) + "," + std::to_string(gradeMv2) + "," + std::to_string(gradeMv3) + ");\n");
    else
        return std::string("\t\trightContractionProductHomogeneous<T>(mv1,mv2,mv3," + std::to_string(gradeMv1) + "," + std::to_string(gradeMv2) + "," + std::to_string(gradeMv3) + ");\n");
}

std::string generateInnerRecursiveFloat64(const unsigned int gradeMv1, const unsigned int gradeMv2, const unsigned int gradeMv3) {

    // call the fine recursive product (right or left contraction)
    if(gradeMv3 == (gradeMv2-gradeMv1))
        return std::string("\t\tleftContractionProductHomogeneous<double>(mv1,mv2,mv3," + std::to_string(gradeMv1) + "," + std::to_string(gradeMv2) + "," + std::to_string(gradeMv3) + ");\n");
    else
        return std::string("\t\trightContractionProductHomogeneous<double>(mv1,mv2,mv3," + std::to_string(gradeMv1) + "," + std::to_string(gradeMv2) + "," + std::to_string(gradeMv3) + ");\n");
}


//// return a string containing all the explicit per-grade inner product functions,case where we handle EUCLIDEAN space (without any basis transformations)
//// the metric has to be diagonal, its coefficient are not all necessarely 1s.
std::string generateInnerExplicit_cpp(const MetaData &metaData,
                                      const ProductTools& product){

    // string that will contain the result
    std::string outputString = "";

    // loop over all multivectors whose grade range from 0 to the dimension
    for(unsigned int gradeMv1=0; gradeMv1<=metaData.dimension; ++gradeMv1){

        // loop over all multivectors mv2 whose grade range from 0 to the dimension
        for(unsigned int gradeMv2=0; gradeMv2<=metaData.dimension; ++gradeMv2){

            // final grade
            auto gradeMv3 = (unsigned int)std::abs((float)gradeMv1-gradeMv2);

            // generate the comments of the current function
            outputString += innerProductExplicitComments(gradeMv1, gradeMv2, gradeMv3);

            // generate the prototype of the current function
            outputString += "template<typename T>\n\tvoid inner_"+std::to_string(gradeMv1)+"_"+std::to_string(gradeMv2)+"(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){\n";;

            // if Mv1 or Mv2 is a scalar, then use the eigen vector to scalar product
            if( (gradeMv1==0) ^ (gradeMv2==0) ){
                if(gradeMv1==0) outputString += "\t\tmv3 += mv1.coeff(0)*mv2;\n";
                else outputString += "\t\tmv3 += mv1*mv2.coeff(0);\n";

            } else { // not a scalar to mv product. Do a regular optimization

                // if one of the k-vector mv1 or mv2 is higher than the threshold, we use the recursive function.
                if( (bin_coeff(metaData.dimension,gradeMv1) > metaData.maxDimPrecomputedProducts) || (bin_coeff(metaData.dimension,gradeMv2) > metaData.maxDimPrecomputedProducts) ){

                    // generate a call to the fine recursive outer product
                    if(metaData.inputMetricDiagonal)  // without basis change
                        outputString += generateInnerRecursive(gradeMv1, gradeMv2, gradeMv3);
                    else  // with basis change
                        outputString += generateInnerRecursiveBasisChange(gradeMv1, gradeMv2, gradeMv3, metaData.dimension);
                }else {

                    // generate all products for this 'gradeMv3'
                    std::list<productComponent<double>> listOfProducts = product.generateExplicitInnerProductListEuclideanSpace(
                            gradeMv1,
                            gradeMv2,
                            metaData.diagonalMetric);


                    // convert all the computed products as a string containing c++ instructions
                    outputString += productListToString(listOfProducts);
                }
            }

            // end of the current function, finishing with a "}"
            outputString += "\t}\n\n\n\t";
        }
    }

    return outputString;
}





//// return a string containing all the explicit per-grade inner product functions, 
std::string generateInnerExplicitBasisChange_cpp(const MetaData &metaData,
                                                 std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor> > transformationMatrices,
                                                 std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor> > inverseTransformationMatrices,
                                                 const ProductTools& product){

    // string that will contain the result
    std::string outputString = "";

    // Outer product 




    // loop over all multivectors whose grade range from 0 to the dimension
    for(unsigned int gradeMv1=0; gradeMv1<=metaData.dimension; ++gradeMv1){

        // loop over all multivectors mv2 whose grade range from 0 to the dimension
        for(unsigned int gradeMv2=0; gradeMv2<=metaData.dimension; ++gradeMv2){

            // final grade
            auto gradeMv3 = (unsigned int)std::abs((float)gradeMv1-gradeMv2);

            // generate the comments of the current function
            outputString += innerProductExplicitComments(gradeMv1, gradeMv2, gradeMv3);

            // generate the prototype of the current function
            outputString += "template<typename T>\n\tvoid inner_"+std::to_string(gradeMv1)+"_"+std::to_string(gradeMv2)+"(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){\n";;

            // if Mv1 or Mv2 is a scalar, then use the eigen vector to scalar product
            if( (gradeMv1==0) ^ (gradeMv2==0) ){
                if(gradeMv1==0) outputString += "\t\tmv3 += mv1.coeff(0)*mv2;\n";
                else outputString += "\t\tmv3 += mv1*mv2.coeff(0);\n";

            } else { // not a scalar to mv product. Do a regular optimization

                // if one of the k-vector mv1 or mv2 is higher than the threshold, we use the recursive function.
                if( (bin_coeff(metaData.dimension,gradeMv1) > metaData.maxDimPrecomputedProducts) || (bin_coeff(metaData.dimension,gradeMv2) > metaData.maxDimPrecomputedProducts) ){

                    // generate a call to the fine recursive outer product
                    if(metaData.inputMetricDiagonal)  // without basis change
                        outputString += generateInnerRecursive(gradeMv1, gradeMv2, gradeMv3);
                    else  // with basis change
                        outputString += generateInnerRecursiveBasisChange(gradeMv1, gradeMv2, gradeMv3, metaData.dimension);
                }else {

                    // generate all products for this 'gradeMv3'
                    std::list<productComponent<double>> listOfProducts = product.generateExplicitInnerProductList(
                            gradeMv1,
                            gradeMv2,
                            inverseTransformationMatrices[gradeMv1],
                            inverseTransformationMatrices[gradeMv2],
                            transformationMatrices[gradeMv3],
                            metaData.diagonalMetric);


                    // convert all the computed products as a string containing c++ instructions
                    outputString += productListToString(listOfProducts);
                }
            }

            // end of the current function, finishing with a "}"
            outputString += "\t}\n\n\n\t";
        }
    }

    return outputString;
}



// generate an array of function pointers referring to each explicit products.
// Start with the declaration of the 2D array of functions pointer, then fill it.
// Note that the size of this functions pointer is (dimension+1)*(dimension+1) because a mv can have (dimension+1) kvectors
std::string generateInnerExplicitFunctionsPointer(const unsigned int &dimension){

    // string that will contain the result
    std::string outputString = "template<typename T>\n"
                                       "\tstd::array<std::array<std::function<void(const Eigen::Matrix<T, Eigen::Dynamic, 1> & , "
                                       "const Eigen::Matrix<T, Eigen::Dynamic, 1> & , Eigen::Matrix<T, Eigen::Dynamic, 1>&)>, "+ std::to_string(dimension+1)+">, "+ std::to_string(dimension+1)+"> innerFunctionsContainer = {{\n";

    // create an entry for each couple of (gradeMv1,gradeMv2)
    for(unsigned int gradeMv1=0; gradeMv1<=dimension; ++gradeMv1){

        // beginning of a line
        outputString += "\t\t{{";

        for(unsigned int gradeMv2=0; gradeMv2<=dimension; ++gradeMv2){

            // create the entry for the couple of (gradeMv1,gradeMv2)
            outputString += "inner_"+std::to_string(gradeMv1)+"_"+std::to_string(gradeMv2)+"<T>";

            // Close the array when gradeMv2 reaches 'dimension'
            // else if gradeMv1 is 'dimension', then put a carriage return, else put a comma
            outputString += (gradeMv2<dimension)?",":"}}"+((gradeMv1==dimension)?std::string("\n"):std::string(",\n"));
        }
    }

    // when all the pointers are initialized, close the array of array with a double }
    outputString += "\t}};";

    return outputString;
}






/*******************************************************************************************************
******************************** EXPLICIT GEOMETRIC PRODUCT DEFINITION ************************************
*******************************************************************************************************/

std::string geometricProductExplicitComments(const unsigned int gradeMv1, const unsigned int gradeMv2, const unsigned int gradeMvResult){
    std::string outerFunctionComments = "/// \\brief Compute the geometric product between two homogeneous multivectors mv1 (grade "+ std::to_string(gradeMv1)+") and mv2 (grade " + std::to_string(gradeMv2) + "). \n\t";
    outerFunctionComments += "/// \\tparam the type of value that we manipulate, either float or double or something else.\n\t";
    outerFunctionComments += "/// \\param mv1 - the first homogeneous multivector of grade " +  std::to_string(gradeMv1)  + " represented as an Eigen::VectorXd\n\t";
    outerFunctionComments += "/// \\param mv2 - the second homogeneous multivector of grade " + std::to_string(gradeMv2) + " represented as a Eigen::VectorXd\n\t";
    outerFunctionComments += "/// \\param mv3 - the result of mv1 mv2 whose grade is "+std::to_string(gradeMvResult)+"\n\t";
    return outerFunctionComments;
}


std::string generateGeometricRecursive(const unsigned int gradeMv1, const unsigned int gradeMv2, const unsigned int gradeMv3) {

    // call the fine recursive product (geometric product)
    return std::string("\t\tgeoProduct<T>(mv1,mv2,mv3," + std::to_string(gradeMv1) + "," + std::to_string(gradeMv2) + "," + std::to_string(gradeMv3) + ");\n");
}

std::string generateGeometricRecursiveFloat64(const unsigned int gradeMv1, const unsigned int gradeMv2, const unsigned int gradeMv3) {

    // call the fine recursive product (geometric product)
    return std::string("\t\tgeoProduct<double>(mv1,mv2,mv3," + std::to_string(gradeMv1) + "," + std::to_string(gradeMv2) + "," + std::to_string(gradeMv3) + ");\n");
}



//// return a string containing all the explicit per-grade inner product functions
std::string generateGeometricExplicit_cpp(const MetaData &metaData,
                                          const ProductTools& product){

    // string that will contain the result
    std::string outputString = "";

    // for all pair of grades
    for(unsigned int gradeMv1=0; gradeMv1<=metaData.dimension; ++gradeMv1){
        for(unsigned int gradeMv2=0; gradeMv2<=metaData.dimension; ++gradeMv2){

            // do not create explicit product for too big multivectors (then recursive calls)
            if( (bin_coeff(metaData.dimension,gradeMv1)>metaData.maxDimPrecomputedProducts) || (bin_coeff(metaData.dimension,gradeMv2)>metaData.maxDimPrecomputedProducts) ){
                    // generate a call to the fine recursive outer product
                    // if(metaData.inputMetricDiagonal)  // without basis change
                    //     outputString += generateGeometricRecursive(gradeMv1, gradeMv2, gradeMv3);
                    // else  // with basis change
                    //     outputString += generateGeometricRecursiveBasisChange(gradeMv1, gradeMv2, gradeMv3, metaData.dimension);

                continue;
            }

            // generate all products for this 'resultingGrade'
            std::vector<std::list<productComponent<double>>> listOfProducts = product.generateExplicitGeometricProductListEuclideanSpace(gradeMv1,
                                                                                                                                         gradeMv2,
                                                                                                                                         metaData.diagonalMetric);

            unsigned int dimension = metaData.dimension;

            // the resulting grades are (possibly) not unique.
            // for each homogeneous multivectors whose grade is neither the (gradeMv1+gradeMv2) nor |gradeMv1-gradeMv2|
            for(unsigned int gradeMv3 = 0;gradeMv3<dimension;++gradeMv3){

                // process only the required products
                if(!(listOfProducts[gradeMv3].empty())){

                    // generate the comments of the current function
                    outputString += geometricProductExplicitComments(gradeMv1, gradeMv2,gradeMv3);

                    // generate the prototype of the current function
                    outputString += "template<typename T>\n\tvoid geometric_"+std::to_string(gradeMv1)+"_"+std::to_string(gradeMv2)+"_"+std::to_string(gradeMv3)+"(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){\n";;

                    // convert all the computed products as a string containing c++ instructions
                    outputString += productListToString(listOfProducts[gradeMv3]);

                    // end of the current function, finishing with a "}"
                    outputString += "\t}\n\n\n\t";
                }
            }
        }
    }

    return outputString;
}


//// return a string containing all the explicit per-grade geometric product functions
std::string generateGeometricExplicitBasisChange_cpp(const MetaData &metaData,
                                                     std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor>> transformationMatrices,
                                                     std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor>> inverseTransformationMatrices,
                                                     const ProductTools& product){

    // string that will contain the result
    std::string outputString = "";

    // for all pair of grades
    for(unsigned int gradeMv1=0; gradeMv1<=metaData.dimension; ++gradeMv1){
        for(unsigned int gradeMv2=0; gradeMv2<=metaData.dimension; ++gradeMv2){

            // do not create explicit product for too big multivectors (then recursive calls)
            if( (bin_coeff(metaData.dimension,gradeMv1)>metaData.maxDimPrecomputedProducts) || (bin_coeff(metaData.dimension,gradeMv2)>metaData.maxDimPrecomputedProducts) )
                continue;

            // generate all products for this 'resultingGrade'
            std::vector<std::list<productComponent<double>>> listOfProducts = product.generateExplicitGeometricProductList(gradeMv1,
                                                                                                                           gradeMv2,
                                                                                                                           inverseTransformationMatrices[gradeMv1],
                                                                                                                           inverseTransformationMatrices[gradeMv2],
                                                                                                                           transformationMatrices,
                                                                                                                           metaData.diagonalMetric);

            unsigned int dimension = transformationMatrices[1].cols();

            // the resulting grades are (possibly) not unique.
            // for each homogeneous multivectors whose grade is neither the (gradeMv1+gradeMv2) nor |gradeMv1-gradeMv2|
            for(unsigned int gradeMv3 = 0;gradeMv3<dimension;++gradeMv3){

                // process only the required products
                if(!(listOfProducts[gradeMv3].empty())){

                    // generate the comments of the current function
                    outputString += geometricProductExplicitComments(gradeMv1, gradeMv2,gradeMv3);

                    // generate the prototype of the current function
                    outputString += "template<typename T>\n\tvoid geometric_"+std::to_string(gradeMv1)+"_"+std::to_string(gradeMv2)+"_"+std::to_string(gradeMv3)+"(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){\n";;

                    // convert all the computed products as a string containing c++ instructions
                    outputString += productListToString(listOfProducts[gradeMv3]);

                    // end of the current function, finishing with a "}"
                    outputString += "\t}\n\n\n\t";
                }
            }
        }
    }

    return outputString;
}




// generate an array of function pointers referring to each explicit products.
// Start with the declaration of the 3D array of functions pointer, then fill it.
// Note that the size of this functions pointer is (dimension+1)*(dimension+1)*(dimension+1)
std::string generateGeometricExplicitFunctionsPointer(const MetaData &metaData){

    // string that will contain the result
    std::string outputString = "template<typename T>\n"
                                       "\tstd::array<std::array<std::array<std::function<void(const Eigen::Matrix<T, Eigen::Dynamic, 1> & , "
                                       "const Eigen::Matrix<T, Eigen::Dynamic, 1> & , Eigen::Matrix<T, Eigen::Dynamic, 1>&)>, "
                               + std::to_string(metaData.dimension+1) + ">, " + std::to_string(metaData.dimension+1) + ">, "
                               + std::to_string(metaData.dimension+1) + "> geometricFunctionsContainer =\n\t{{\n";

    // for each triplet of grade
    for(unsigned int gradeMv3=0; gradeMv3<=metaData.dimension; ++gradeMv3) {

        // start an function pointer array for Mv3 (even if it is empty)
        outputString += "\t\t{{\n";

        for(unsigned int gradeMv1=0; gradeMv1<=metaData.dimension; ++gradeMv1) {

            // start an function pointer array for Mv2 (even if it is empty)
            outputString += "\t\t\t{{";

            for(unsigned int gradeMv2=0; gradeMv2<=metaData.dimension; ++gradeMv2) {

                // compute outer / inner / geometric grades
                unsigned int gradeOuter = gradeMv1 + gradeMv2;
                unsigned int gradeInner = (unsigned int) std::abs(int(gradeMv1)-int(gradeMv2));
                unsigned int gradeGeometric = (gradeMv3 - gradeInner)%2;
                unsigned int gradeMax = (2*metaData.dimension)-gradeOuter;

                // ignore the outer and inner products
                // the grade of the result has to be in [gradeMv1-gradeMv2, gradeMv1-gradeMv2+2, gradeMv1-gradeMv2+4,...]
                if((gradeMv3 > gradeInner) && (gradeMv3 <= gradeMax) && (gradeGeometric == 0) && (gradeMv3 < gradeOuter) && ( bin_coeff(metaData.dimension,gradeMv1) <= metaData.maxDimPrecomputedProducts) && (bin_coeff(metaData.dimension,gradeMv2) <= metaData.maxDimPrecomputedProducts) ){
                    // create the entry for the triple of (gradeMv3,gradeMv1,gradeMv2)
                    outputString += "geometric_" + std::to_string(gradeMv1) + "_" + std::to_string(gradeMv2)+ "_" + std::to_string(gradeMv3) + "<T>";
                }else{
                    outputString += "{}";
                }

                // Close the array when gradeMv2 reaches 'dimension'
                // else if gradeMv1 is 'dimension', then put a carriage return, else put a comma
                outputString += (gradeMv2 < metaData.dimension) ? "," : "}}" + ((gradeMv1 == metaData.dimension) ? std::string("\n") : std::string(",\n"));
            }
        }
        outputString += (gradeMv3 < metaData.dimension) ? "\t\t}},\n" : "\t\t}}\n";
    }

    // when all the pointers are initialized, close the array of array with a double }
    outputString += "\t}};";

    return outputString;
}


std::string recursiveGeometricProductCallToString(const unsigned int maxSize, bool changeMetricToOrtrhogonal){
    std::string output = "\n                // choose iterative or recursive\n";
    output += "                if( (itMv1.vec.size()>=" + std::to_string(maxSize)
              + ") || (itMv2.vec.size()>=" + std::to_string(maxSize)
              + ") ){\n";

    if(changeMetricToOrtrhogonal){
        output += "                   Mvec<T> tmp;\n";
        output += "                   // geometric product from mv1 and mv2 in the orthogonal basis\n";
        output += "                   geoProduct<T>(transformationMatricesInverse<T>[itMv1.grade]*itMv1.vec, transformationMatricesInverse<T>[itMv2.grade]*itMv2.vec, tmp, itMv1.grade, itMv2.grade, itMv1.grade + itMv2.grade);\n";
        output += "                   for(auto & itMvTmp : tmp.mvData){\n";
        output += "                       auto it = mv3.createVectorXdIfDoesNotExist(itMvTmp.grade);\n";
        output += "                       // back transform the result in the original basis\n";
        output += "                       it->vec += transformationMatrices<T>[itMvTmp.grade] * itMvTmp.vec;\n";
        output += "                       if(!((it->vec.array() != 0.0).any())){\n";
        output += "                           mv3.mvData.erase(it);\n                       }\n                   }\n";

    }else{
        output += "                   geoProduct<T>(itMv1.vec, itMv2.vec, mv3, itMv1.grade, itMv2.grade, itMv1.grade + itMv2.grade);\n";
    }

    output += "                   continue;\n                }\n";

    return output;
}


std::string recursiveGeometricProductCallToStringFloat64(const unsigned int maxSize, bool changeMetricToOrtrhogonal){
    std::string output = "\n                // choose iterative or recursive\n";
    output += "                if( (itMv1.vec.size()>=" + std::to_string(maxSize)
              + ") || (itMv2.vec.size()>=" + std::to_string(maxSize)
              + ") ){\n";

    if(changeMetricToOrtrhogonal){
        output += "                   Mvec<double> tmp;\n";
        output += "                   // geometric product from mv1 and mv2 in the orthogonal basis\n";
        output += "                   geoProduct<double>(transformationMatricesInverse<double>[itMv1.grade]*itMv1.vec, transformationMatricesInverse<double>[itMv2.grade]*itMv2.vec, tmp, itMv1.grade, itMv2.grade, itMv1.grade + itMv2.grade);\n";
        output += "                   for(auto & itMvTmp : tmp.mvData){\n";
        output += "                       auto it = mv3.createVectorXdIfDoesNotExist(itMvTmp.grade);\n";
        output += "                       // back transform the result in the original basis\n";
        output += "                       it->vec += transformationMatrices<double>[itMvTmp.grade] * itMvTmp.vec;\n";
        output += "                       if(!((it->vec.array() != 0.0).any())){\n";
        output += "                           mv3.mvData.erase(it);\n                       }\n                   }\n";
    }else{
        output += "                   geoProduct<double>(itMv1.vec, itMv2.vec, mv3, itMv1.grade, itMv2.grade, itMv1.grade + itMv2.grade);\n";
    }

    output += "                   continue;\n                }\n";

    return output;
}
















