#include <fstream>
#include <vector>
#pragma once
// This class is used to store the bins for a histogram
template <typename T>
class Bins {
    public:
        std::vector<int> shape;
        std::vector<T> bins;

        void addFeature(int nBins);
        void addCount(std::vector<int> indices, T count);
        T getCount(std::vector<int> indices);
        int getIndex(std::vector<int> indices);
        std::vector<std::vector<int>> getIndices(std::vector<int> shape);

        void clear();

};

// This class is used realize the structure of a histogram and to store the data
template <typename T>
class Histogram {
    public:
        Bins<T> bins;
        std::vector<std::vector<double>> bin_minmax;
        std::vector<std::string> feature_names;

        void addFeature(std::string name, int nBins, double min, double max);
        void addEvent(std::vector<double> features, T count = 1);
        void saveData(std::string filename);

        std::vector<std::vector<int>> getIndices();
        std::vector<double> getArgumentValues(std::vector<int> indices);
        double binValue(int bin, int feature);

        void print();
        void clear();
};

////////////////////////////////////
// HISTOGRAM CLASS IMPLEMENTATION //
////////////////////////////////////

template <typename T>
void Histogram<T>::addFeature(std::string name, int nBins, double min, double max) {
    bins.addFeature(nBins);
    bin_minmax.push_back({min, max});
    feature_names.push_back(name);
}
template <typename T>
void Histogram<T>::addEvent(std::vector<double> features, T count) {
    std::vector<int> indices;
    for (int i = 0; i < features.size(); i++) {
        int bin = (features[i] - bin_minmax[i][0]) / (bin_minmax[i][1] - bin_minmax[i][0]) * bins.shape[i];
        if (bin == bins.shape[i]) {
            bin -= 1;
        }
        indices.push_back(bin);
    }
    bins.addCount(indices, count);
}
template <typename T>
void Histogram<T>::saveData(std::string filename) {
    std::ofstream file;
    file.open(filename);
    // Write the header
    for (int i = 0; i < feature_names.size(); i++) {
        file << feature_names[i] << ",";
    }
    file << "count\n";
    // Write the data
    std::vector<std::vector<int>> indexList = bins.getIndices(bins.shape);
    for (int i = 0; i < indexList.size(); i++) {
        for (int j = 0; j < indexList[i].size(); j++) {
            file << binValue(indexList[i][j], j) << ",";
        }
        file << bins.getCount(indexList[i]) << "\n";
    }
    file.close();
}
template <typename T>
double Histogram<T>::binValue(int binIndex, int featureIndex) {
    double min = bin_minmax[featureIndex][0];
    double max = bin_minmax[featureIndex][1];
    double bin_size = (max - min) / bins.shape[featureIndex];
    return min + bin_size * (binIndex + 0.5);
}
template <typename T>
std::vector<std::vector<int>> Histogram<T>::getIndices() {
    return bins.getIndices(bins.shape);
}
template <typename T>
std::vector<double> Histogram<T>::getArgumentValues(std::vector<int> indices) {
    std::vector<double> values;
    int bin = bins.getIndex(indices);
    for (int feature = 0; feature < indices.size(); feature++) {
        values.push_back(binValue(bin, feature));
    }
    return values;
}
template <typename T>
void Histogram<T>::print() {
    std::vector<std::vector<int>> indexList = bins.getIndices(bins.shape);
    for (int i = 0; i < indexList.size(); i++) {
        for (int j = 0; j < indexList[i].size(); j++) {
            std::cout << binValue(indexList[i][j], j) << ",";
        }
        std::cout << bins.getCount(indexList[i]) << "\n";
    }
}

template <typename T>
void Histogram<T>::clear() {
    bins.clear();
}

////////////////////////////////////
//// BINS CLASS IMPLEMENTATION /////
////////////////////////////////////

template <typename T>
void Bins<T>::addFeature(int nBins) {
    if (shape.size() == 0) {
        shape.push_back(nBins);
        bins = std::vector<T>(nBins);
    } else {
        shape.push_back(nBins);
        int old_size = bins.size();
        for (int i = 0; i < old_size * nBins - old_size; i++) {
            bins.push_back(0);
        }
    }
}
template <typename T>
void Bins<T>::addCount(std::vector<int> indices, T count) {
    int index = getIndex(indices);
    bins[index] += count;
}
template <typename T>
T Bins<T>::getCount(std::vector<int> indices) {
    int index = getIndex(indices);
    return bins[index];
}
template <typename T>
int Bins<T>::getIndex(std::vector<int> indices) {
    int index = 0;
    int buffer = bins.size();
    for (int i = 0; i < indices.size(); i++) {
        buffer /= shape[i];
        index += indices[i] * buffer;
    }
    return index;
}
template <typename T>
std::vector<std::vector<int>> Bins<T>::getIndices(std::vector<int> shape) {
    if (shape.size() == 1) {
        std::vector<std::vector<int>> indexList;
        for (int i = 0; i < shape[0]; i++) {
            indexList.push_back({i});
        }
        return indexList;
    } else {
        std::vector<std::vector<int>> indexList;
        std::vector<int> sub_shape(shape.begin() + 1, shape.end());
        std::vector<std::vector<int>> sub_indices = getIndices(sub_shape);
        for (int i = 0; i < shape[0]; i++) {
            for (int j = 0; j < sub_indices.size(); j++) {
                std::vector<int> index = {i};
                index.insert(index.end(), sub_indices[j].begin(), sub_indices[j].end());
                indexList.push_back(index);
            }
        }
        return indexList;
    }
}

template <typename T>
void Bins<T>::clear() {
    bins.clear();
    shape.clear();
}