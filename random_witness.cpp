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

    // set up the output file
    std::ofstream outfile;
    outfile.open("random_witness.csv");
    outfile << "I_trMrho_e,I_sign_e\n";

    // simulation parameters
    // int nFunctionals = 10000;
    // int nSamples = 20000;
    int nFunctionals = 100;
    int nSamples = 1000;

    // variables
    double trMrho;
    double entanglement;
    double I_tr_e;
    double I_sign_e;
    Eigen::Matrix4cd M;
    Eigen::Matrix4cd rho;
    JointProbability joint;
    
    for (int f = 0; f < nFunctionals; f++) {
        M = Linalg::random4Witness();
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
        I_tr_e = joint.MutualInformation(0, 1);
        I_sign_e = joint.WitnessInformation(0, 1);
        outfile << I_tr_e << "," << I_sign_e << "\n";

        std::cout << "nSamples = " << nSamples << std::endl;
        std::cout << "I(tr(Mrho) : entanglement) = " << I_tr_e << std::endl;
        std::cout << "I(sign : entanglement) = " << I_sign_e << "\n\n";

        hist.clear();
        joint.clear();
    }
    outfile.close();
    return 0;
}