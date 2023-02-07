#include <iostream>
#include <fstream>

#include "matrix_loader.h"
#include "histogram_tools.h"
#include "information_tools.h"
#include "linalg_tools.h"



int main (int argc, char *argv[]) {

    Histogram<int> hist;
    Histogram<int> MI;
    MI.addFeature("I_trMrho_e", 100, 0, 1);
    MI.addFeature("I_sign_e", 100, 0, 1);

    int nFunctionals = 1;
    int nSamples = 100000;
    double trMrho;
    double entanglement;
    double I_tr_e;
    double I_sign_e;
    Eigen::Matrix4cd M;
    Eigen::Matrix4cd rho;
    JointProbability joint;
    
    for (int f = 0; f < nFunctionals; f++) {
        hist = Histogram<int>();
        hist.addFeature("tr(Mrho)", 1000, -1, 1);
        hist.addFeature("entanglement", 2, -0.5, 1.5);
        M = Linalg::random4State();
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
        MI.addEvent({I_tr_e, I_sign_e});
        std::cout << "\nI(tr(Mrho) : entanglement) = " << I_tr_e << std::endl;
        std::cout << "I(sign : entanglement) = " << I_sign_e << std::endl;
    }
    

    return 0;
}