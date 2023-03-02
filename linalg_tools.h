#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <complex>

std::vector<Eigen::Matrix<std::complex<double>, 16,4>> partialTraceMatricesRight {
    Eigen::Matrix<std::complex<double>, 16, 4>{
    {1,0,0,0},
    {0,0,0,0},
    {0,0,0,0},
    {0,0,0,0},
    {0,1,0,0},
    {0,0,0,0},
    {0,0,0,0},
    {0,0,0,0},
    {0,0,1,0},
    {0,0,0,0},
    {0,0,0,0},
    {0,0,0,0},
    {0,0,0,1},
    {0,0,0,0},
    {0,0,0,0},
    {0,0,0,0},
    },
    Eigen::Matrix<std::complex<double>, 16, 4>{
    {0,0,0,0},
    {1,0,0,0},
    {0,0,0,0},
    {0,0,0,0},
    {0,0,0,0},
    {0,1,0,0},
    {0,0,0,0},
    {0,0,0,0},
    {0,0,0,0},
    {0,0,1,0},
    {0,0,0,0},
    {0,0,0,0},
    {0,0,0,0},
    {0,0,0,1},
    {0,0,0,0},
    {0,0,0,0},
    },
    Eigen::Matrix<std::complex<double>, 16, 4>{
    {0,0,0,0},
    {0,0,0,0},
    {1,0,0,0},
    {0,0,0,0},
    {0,0,0,0},
    {0,0,0,0},
    {0,1,0,0},
    {0,0,0,0},
    {0,0,0,0},
    {0,0,0,0},
    {0,0,1,0},
    {0,0,0,0},
    {0,0,0,0},
    {0,0,0,0},
    {0,0,0,1},
    {0,0,0,0},
    },
    Eigen::Matrix<std::complex<double>, 16, 4>{
    {0,0,0,0},
    {0,0,0,0},
    {0,0,0,0},
    {1,0,0,0},
    {0,0,0,0},
    {0,0,0,0},
    {0,0,0,0},
    {0,1,0,0},
    {0,0,0,0},
    {0,0,0,0},
    {0,0,0,0},
    {0,0,1,0},
    {0,0,0,0},
    {0,0,0,0},
    {0,0,0,0},
    {0,0,0,1},
    },
};

std::vector<Eigen::Matrix<std::complex<double>, 4, 16>> partialTraceMatricesLeft {
    Eigen::Matrix<std::complex<double>, 4, 16>{
    {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0},
    },
    Eigen::Matrix<std::complex<double>, 4, 16>{
    {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0},
    },
    Eigen::Matrix<std::complex<double>, 4, 16>{
    {0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0},
    },
    Eigen::Matrix<std::complex<double>, 4, 16>{
    {0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    }
};

namespace Linalg{

        std::complex<double> randomComplex() {
            double subNormalization = rand() / (double)RAND_MAX;
            double realPart = rand() / (double)RAND_MAX;
            std::complex<double> randomNumber = std::complex<double>(sqrt(realPart) , sqrt(1 - realPart));
            return randomNumber * subNormalization;
        }
        void swapIJ (Eigen::Matrix4cd& rho, int i, int j, int k, int l) {
            // Swap the elements (i,j) and (k,l) of the matrix rho
            std::complex<float> temp;
            temp = rho(i, j);
            rho(i, j) = rho(k, l);
            rho(k, l) = temp;
        }
        Eigen::Matrix4cd partialTranpose2nd (Eigen::Matrix4cd rho) {
            // Partial transpose the second qubit
            Eigen::Matrix4cd rho2 = rho;
            swapIJ(rho2, 0, 1, 1, 0);
            swapIJ(rho2, 2, 1, 3, 0);
            swapIJ(rho2, 1, 2, 0, 3);
            swapIJ(rho2, 3, 2, 2, 3);
            return rho2;
        }
        bool isEntangledPPT(Eigen::Matrix4cd rho) {
            // Check if the state is entangled through the PPT criterion
            bool result = partialTranpose2nd(rho).eigenvalues().real().minCoeff() < 0;
            return result;
        }
        Eigen::MatrixXcd random4State() {
            // Start with an empty 16x1 matrix
            Eigen::MatrixXcd psi{{0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0}};
            std::complex<double> norm = 0;
            
            // Fill it with random complex numbers
            for(int i = 0; i < 16; i++) {
                psi(i,0) = randomComplex();
                norm += psi(i,0) * std::conj(psi(i,0));
            }
            // Normalize the state
            for(int i = 0; i < 16; i++) {
                psi(i,0) /= std::sqrt(norm);
            }
            // Calculate the density matrix of the state
            Eigen::MatrixXcd rho{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
            Eigen::MatrixXcd psipsi = psi * psi.adjoint();
            // Return the partial trace of the density matrix over the second qubit
            for(int i = 0; i < 4; i++) {
                rho += partialTraceMatricesLeft[i] * psipsi * partialTraceMatricesRight[i];
            }
            return rho;
        }

        Eigen::Matrix4cd random4Functional() {
            // Generate a random 4x4 hermitian matrix that has norm = 1
            Eigen::Matrix4cd M = Eigen::Matrix4cd::Random();
            Eigen::Matrix4cd N = M.adjoint();
            M = M + N;
            M = M / sqrt(((M * M).trace()));
            return M;
        }
        Eigen::Matrix4cd random4Witness() {
            std::complex<double> alpha, beta, gamma;

            alpha = randomComplex();
            beta = randomComplex();
            gamma = std::complex<double>(rand() / (double)RAND_MAX, 0);

            if (rand() % 2 == 0) {
                gamma = -gamma;
            }

            std::complex<double> one(1, 0);
            
            // Make sure the witness has norm = 1
            Eigen::Matrix4cd M;
            M(0, 0) = one + gamma;
            M(1, 1) = one - gamma;
            M(2, 2) = one - gamma;
            M(3, 3) = one + gamma;

            M(0, 3) = alpha + beta;
            M(3, 0) = std::conj(M(0, 3));
            M(1, 2) = alpha - beta;
            M(2, 1) = std::conj(M(1,2));

            M = M / (M * M).trace();
            return M;
        }

};
