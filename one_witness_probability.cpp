#include <iostream>
#include <fstream>

#include "matrix_loader.h"
#include "histogram_tools.h"
#include "information_tools.h"
#include "linalg_tools.h"

double InnerProductMatrixRho(Eigen::Matrix4cd M, Eigen::Matrix4cd rho) {
    std::complex<double> result = 0;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            result += std::conj(M(i,j)) * rho(i, j);
        }
    }
    return result.real();
}


int main (int argc, char *argv[]) {

    // initialize a histogram
    Histogram<int> hist;

    // simulation parameters
    int nSamples = 20000;

    // variables
    double trMrho;
    double entanglement;
    Eigen::Matrix4cd M;
    Eigen::Matrix4cd rho;
    JointProbability joint;
    
    M = Linalg::random4Witness();

    // Eigen::VectorXcd eigenvals = M.eigenvalues();

    // std::cout << "eigenvals = \n" << eigenvals << std::endl;

    // WORKING WITNESS
    // M(0, 0) = std::complex<double>(0.5, 0);
    // // M(1, 1) = 1;
    // // M(2, 2) = 1;
    // M(3, 3) = std::complex<double>(0.5, 0);
    // M(1, 2) = std::complex<double>(0.5, 0);
    // M(2, 1) = std::complex<double>(0.5, 0);
    // M = M / sqrt(((M * M).trace()));

    hist = Histogram<int>();
    hist.addFeature("tr(Mrho)", 100, -1, 1);
    hist.addFeature("entanglement", 2, -0.5, 1.5);

    for (int i = 0; i < nSamples; i++) {
        rho = Linalg::random4State();
        trMrho = InnerProductMatrixRho(M, rho);
        entanglement = Linalg::isEntangledPPT(rho);
        if (entanglement) {
            hist.addEvent({trMrho, 1});
        } else {
            hist.addEvent({trMrho, 0});
        }
    }

    joint = JointProbability(hist);
    hist.saveData("one_witness_histogram.csv");
    joint.saveData("one_witness_distribution.csv");

    hist.clear();
    joint.clear();

    return 0;
}
