#include "histogram_tools.h"
#include <vector>


class JointProbability : public Histogram<double>{
    // This stores the joint probability of the histogram provided in the constructor
    // without regarding the values of each outcome
    public:
        JointProbability() {};
        JointProbability(Histogram<int> histogram);
        JointProbability(std::vector<int> shape) {
            for (int nBins : shape) {
                bins.addFeature(nBins);
            }
        };
        
        double P(std::vector<int> indices);
        double PX(int X, int x);
        double MutualInformation(int X, int Y);
        double WitnessInformation (int X, int Y, double zero);
};

JointProbability::JointProbability(Histogram<int> histogram) {
    std::vector<int> shape = histogram.bins.shape;
    for (int nBins : shape) {
        bins.addFeature(nBins);
    }
    bin_minmax = histogram.bin_minmax;
    feature_names = histogram.feature_names;

    int total_counts = 0;
    for (int counts : histogram.bins.bins) {
        total_counts += counts;
    }

    std::vector<std::vector<int>> indexList = getIndices();
    for (int i = 0; i < indexList.size(); i++) {
        bins.addCount(indexList[i], (double) histogram.bins.getCount(indexList[i]) / total_counts);
    }
}
double JointProbability::P(std::vector<int> indices) {
    return bins.getCount(indices);
}

double JointProbability::PX(int X, int x) {
    // This calculates the probability of a single feature
    // given the joint probability
    // P(X=x) = sum_y P(X=x, Y=y)
    double sum = 0;
    std::vector<std::vector<int>> indexList = getIndices();
    for (std::vector<int> index : indexList) {
        if (index[X] == x) {
            sum += P(index);
        }
    }
    return sum;
}

double JointProbability::MutualInformation(int X, int Y) {
    // This calculates the mutual information between two features
    // I(X;Y) = sum_x sum_y P(X=x, Y=y) log(P(X=x, Y=y) / P(X=x)P(Y=y))
    double sum = 0;
    std::vector<std::vector<int>> indexList = getIndices();
    for (std::vector<int> index : indexList) {
        double PXY = P(index);
        if (PXY != 0) {
            sum += PXY * log2(PXY / (PX(X, index[X]) * PX(Y, index[Y])));
        }
    }
    return sum;
}

double JointProbability::WitnessInformation (int X, int Y, double zero = 0) {
    JointProbability joint_coarse({2, 2});
    joint_coarse.feature_names = {"sign", "entanglement"};
    joint_coarse.bin_minmax = {{-2, 2}, {-0.5, 1.5}};    
    
    std::vector<std::vector<int>> indexList = getIndices();

    for (int i = 0; i < indexList.size(); i++) {
        if (binValue(indexList[i][X], X) > zero) {
            if (binValue(indexList[i][Y], Y) == 1) {
                joint_coarse.bins.addCount({1, 1}, P(indexList[i]));
            } else {
                joint_coarse.bins.addCount({1, 0}, P(indexList[i]));
            }
        } else {
            if (binValue(indexList[i][Y], Y) == 1) {
                joint_coarse.bins.addCount({0, 1}, P(indexList[i]));
            } else {
                joint_coarse.bins.addCount({0, 0}, P(indexList[i]));
            }
        }
    }

    // print();
    // std::cout << "\n";
    // joint_coarse.print();
    // std::cout << joint_coarse.P({0, 0}) << "\n";
    // std::cout << joint_coarse.P({0, 1}) << "\n";
    // std::cout << joint_coarse.P({1, 0}) << "\n";
    // std::cout << joint_coarse.P({1, 1}) << "\n\n";

    // std::cout << joint_coarse.PX(0, 0) << "\n";
    // std::cout << joint_coarse.PX(0, 1) << "\n";
    // std::cout << joint_coarse.PX(1, 0) << "\n";
    // std::cout << joint_coarse.PX(1, 1) << "\n";

    return joint_coarse.MutualInformation(0, 1);
}