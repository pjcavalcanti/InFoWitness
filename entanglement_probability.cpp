#include <iostream>
#include <Eigen/Dense>
#include "linalg_tools.h"
#include "information_tools.h"

int main() {

    // simulation parameters
    int nSamples = 100000;
    int nEntangled = 0;

    Eigen::Matrix4cd rho;
    double entanglement;

    std::vector<Eigen::Matrix4cd> rhos;
    for (int i = 0; i < nSamples; i++) {
        rhos.push_back(Linalg::random4State());
    }

    for (int i = 0; i < nSamples; i++) {
            rho = rhos[i];
            entanglement = Linalg::isEntangledPPT(rho);
        if (entanglement) {
            nEntangled += 1;
        }
    }

    std::cout << (double) nEntangled / nSamples << std::endl;

    return 0;
}