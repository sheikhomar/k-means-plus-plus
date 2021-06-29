#include <clustering/local_search.hpp>

using namespace clustering;

LocalSearch::LocalSearch(uint k, uint s) : numOfClusters(k), swapSize(s)
{
}

std::shared_ptr<ClusteringResult>
LocalSearch::run(const blaze::DynamicMatrix<double> &data)
{
    KMeans kMeansAlg(this->numOfClusters, true, 100U, 0.0001, 42);
    auto result = kMeansAlg.run(data);
    std::cout << "Cluster labels: \n" << result->getCentroids() <<  "\n" ;
    return result;
}
