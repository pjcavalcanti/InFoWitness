#include <iostream>
#include <fstream>

#include "matrix_loader.h"
#include "histogram_tools.h"
#include "information_tools.h"
#include "linalg_tools.h"



int main (int argc, char *argv[]) {

    Histogram<int> hist;
    // Histogram<int> MI;
    // MI.addFeature("I_trMrho_e", 10000, 0, 1);
    // MI.addFeature("I_sign_e", 10000, 0, 1);
    // set up the ouput file
    std::ofstream outfile;
    outfile.open("random_functional.csv");
    outfile << "I_trMrho_e,I_sign_e\n";

    int nFunctionals = 100;
    int nSamples = 20000;
    double trMrho;
    double entanglement;
    double I_tr_e;
    double I_sign_e;
    Eigen::Matrix4cd M;
    Eigen::Matrix4cd rho;
    JointProbability joint;
    
    for (int f = 0; f < nFunctionals; f++) {
        M = Linalg::random4State();
        hist = Histogram<int>();
        hist.addFeature("tr(Mrho)", 100, -1, 1);
        hist.addFeature("entanglement", 2, -0.5, 1.5);
        for (int i = 0; i < nSamples; i++) {
            rho = Linalg::random4State();
            trMrho = (M * rho).trace().real();
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
        std::cout << "nSamples = " << nSamples << std::endl;
        std::cout << "I(tr(Mrho) : entanglement) = " << I_tr_e << std::endl;
        std::cout << "I(sign : entanglement) = " << I_sign_e << "\n\n";
        // MI.addEvent({I_tr_e, I_sign_e});
        outfile << I_tr_e << "," << I_sign_e << "\n";
    }
    // MI.saveData("random_functional.csv");
    outfile.close();
    

    return 0;
}