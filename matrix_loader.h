// This contains a function that loads a 4x4 complex matrix from a file of the form
// (a1,b1) (a2,b2) (a3,b4) (a4,b4)
// (a5,b5) (a6,b6) (a7,b7) (a8,b8)
// (a9,b9) (a10,b10) (a11,b11) (a12,b12)
// (a13,b13) (a14,b14) (a15,b15) (a16,b16)
// where each element is a complex number of the form a+bi represented as (a,b)
#include <Eigen/Dense>
#include <fstream>

Eigen::Matrix4cd load4Matrix (std::string filename){
    // Read the complex matrix M from a file where each line is a row of the matrix
    // and each element is separated by a space, and each element is a complex number of the form a+bi represented as (a,b)
    // The matrix is assumed to be 4x4
    Eigen::Matrix4cd M;
    // Setting up the file
    std::ifstream file;
    std::string line;

    file.open(filename);
    // Setting up buffers
    std::string numberBuffer;
    std::complex <double> Mij;
    // Reading the file
    int row = -1;
    int column = -1;
    while (getline (file, line)) {
        row++;
        numberBuffer = "";
        for (int i = 0; i < line.length(); i++) {
            if (line[i] != ' ' ) {
                numberBuffer += line[i];
            } else {
                column++;
                Mij = std::complex<double>(
                    std::stod(numberBuffer.substr(1, numberBuffer.find(","))),
                    std::stod(numberBuffer.substr(numberBuffer.find(",") + 1, numberBuffer.find(")")))
                    );
                M(row, column) = Mij;
                numberBuffer = "";
            }
        }
        column++;
        Mij = std::complex<double>(
            std::stod(numberBuffer.substr(1, numberBuffer.find(","))),
            std::stod(numberBuffer.substr(numberBuffer.find(",") + 1, numberBuffer.find(")")))
            );
        M(row, column) = Mij;
        column = -1;
    }
    file.close();
    return M;
}