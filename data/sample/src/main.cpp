#include <iostream>
#include <cstdlib>
#include <project_namespace/Mvec.hpp>



int main(){

    // sample instructions
    std::cout << "metric : \n" << project_namespace::metric << std::endl;

    // accessor
    project_namespace::Mvec<double> mv1;
    mv1[project_namespace::scalar] = 1.0;
    mv1[project_namespace::Eproject_first_vector_basis] = 42.0;
    std::cout << "mv1 : " << mv1 << std::endl;

    project_namespace::Mvec<double> mv2;
    mv2[project_namespace::Eproject_first_vector_basis] = 1.0;
    mv2[project_namespace::Eproject_second_vector_basis] = 2.0;
    mv2 += project_namespace::I<double>() + 2*project_namespace::eproject_first_vector_basisproject_second_vector_basis<double>();
    std::cout << "mv2 : " << mv2 << std::endl << std::endl;

    // some products
    std::cout << "outer product     : " << (mv1 ^ mv2) << std::endl;
    std::cout << "inner product     : " << (mv1 | mv2) << std::endl;
    std::cout << "geometric product : " << (mv1 * mv2) << std::endl;
    std::cout << "left contraction  : " << (mv1 < mv2) << std::endl;
    std::cout << "right contraction : " << (mv1 > mv2) << std::endl;
    std::cout << std::endl;

    // some tools
    std::cout << "grade : " << mv1.grade()  << std::endl;
    std::cout << "norm  : " << mv1.norm()  << std::endl;
    mv1.clear();
    if(mv1.isEmpty()) std::cout << "mv1 is empty: ok" << std::endl;

    return EXIT_SUCCESS;
}
