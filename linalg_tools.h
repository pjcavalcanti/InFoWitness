#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <complex>

namespace Linalg{
    // public:
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


        Eigen::Matrix4cd random4State() {
            // Generate a random 4x4 complex positive matrix with trace = 1
            Eigen::Matrix4cd rho = Eigen::Matrix4cd::Random();
            rho = rho * rho.adjoint();
            rho = rho / rho.trace();
            return rho;
        }
        Eigen::Matrix4cd random4Functional() {
            // Generate a random 4x4 complex matrix has abs(trace) = 1
            Eigen::Matrix4cd M = Eigen::Matrix4cd::Random();
            Eigen::Matrix4cd N = M.adjoint();
            M = M + N;
            M = M / sqrt(((M * M).trace()));
            return M;
        }
        Eigen::Matrix4cd random4Witness() {
            std::complex<double> alpha, beta, gamma, norm;
            // alpha = std::complex<double>(rand() / (double)RAND_MAX, rand() / (double)RAND_MAX);
            // beta = std::complex<double>(rand() / (double)RAND_MAX, rand() / (double)RAND_MAX);
            // gamma = std::complex<double>(rand() / (double)RAND_MAX, rand() / (double)RAND_MAX);

            // alpha = alpha / sqrt(alpha * std::conj(alpha)) * std::complex(rand() / (double) RAND_MAX, 0);
            // beta = beta / sqrt(beta * std::conj(beta)) std::complex(rand() / (double) RAND_MAX, 0);
            // gamma = gamma / sqrt(gamma * std::conj(gamma)) std::complex(rand() / (double) RAND_MAX, 0);

            alpha = std::complex<double>(rand() / (double)RAND_MAX, 0);
            beta = std::complex<double>(rand() / (double)RAND_MAX, 0);
            gamma = std::complex<double>(rand() / (double)RAND_MAX, 0);

            if (rand() % 2 == 0) {
                alpha = -alpha;
            }
            if (rand() % 2 == 0) {
                beta = -beta;
            }
            if (rand() % 2 == 0) {
                gamma = -gamma;
            }

            std::complex<double> one(1, 0);
            
            std::cout << "alpha = " << alpha << std::endl;
            std::cout << "beta = " << beta << std::endl;
            std::cout << "gamma = " << gamma << std::endl;
            std::cout << "one = " << one << std::endl;
            // Generate a random 4x4 complex matrix has abs(trace) = 1
            Eigen::Matrix4cd M;
            M(0, 0) = one + gamma;
            M(1, 1) = one - gamma;
            M(2, 2) = one - gamma;
            M(3, 3) = one + gamma;

            M(3, 0) = alpha + beta;
            M(0, 3) = alpha + beta;
            M(2, 1) = alpha - beta;
            M(1, 2) = alpha - beta;

            M = M / (M).trace();
            return M;
        }

};
