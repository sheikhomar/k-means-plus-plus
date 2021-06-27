#include <kmeans/kmeans.hpp>

using namespace kmeans;

KMeans::KMeans(uint k, bool kpp, uint miter, double convDiff) :
    numOfClusters(k),
    initKMeansPlusPlus(kpp),
    maxIterations(miter),
    convergenceDiff(convDiff) {
}

blaze::DynamicVector<size_t> 
KMeans::run(const blaze::DynamicMatrix<double>& data) {
    blaze::DynamicVector<size_t> d = {10, 20, 20};
    return d;
}
