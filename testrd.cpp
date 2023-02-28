#include <iostream>
#include <Eigen/Dense>
#include "linalg_tools.h"

int main() {

    Eigen::MatrixXcd rho(4,4);
    
    for(int i = 0; i < 10; i++) {
        rho = Linalg::random4State();
        std::cout << "rho =" << std::endl;
        std::cout << rho << std::endl;
        std::cout << "trace =" << std::endl;
        std::cout << rho.trace() << std::endl;
        std::cout << "purity =" << std::endl;
        std::cout << (rho * rho).trace() << std::endl;
        std::cout << "eigenvals =" << std::endl;
        std::cout << (rho).eigenvalues() << std::endl;
        std::cout << "\n\n\n" << std::endl; 
    }

    return 0;
}