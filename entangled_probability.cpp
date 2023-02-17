#include <iostream>

#include "linalg_tools.h"

int main(int arc, char*arcv[]) {

    // simulation parameters
    int nSamples = 1000000;
    int nEntangled = 0;

    Eigen::Matrix4cd rho;
    bool entanglement;

    for (int i = 0; i < nSamples; i++) {
        rho = Linalg::random4State();
        entanglement = Linalg::isEntangledPPT(rho);
        if (entanglement) {
            nEntangled += 1;
        }
    }

    std::cout << (double) nEntangled / nSamples;

    return 0;
}
