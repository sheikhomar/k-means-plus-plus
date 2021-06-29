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

    auto clusterAssignments = result->getClusterAssignments();
    auto centers = result->getCentroids();

    auto n = data.rows();
    auto k = clusterAssignments.getNumberOfClusters();

    std::vector<uint> availableIndices(n);
    std::iota(availableIndices.begin(), availableIndices.end(), 0);

    // Find the cost of the clustering done by k-Means.
    double bestCost = clusterAssignments.calcCost();
    auto bestCenters = centers;

    auto swapClusterAssignments = clusterAssignments;

    printf("Cost before swaps %0.5f\n", bestCost);

    std::cout << "Best centers: \n" << bestCenters << "\n";

    for (size_t c = 0; c < k; c++)
    {
        for (size_t p = 0; p < n; p++)
        {
            // Swap one center (c) with a point (p)
            blaze::row(centers, c) = blaze::row(data, p);

            // Reassign points to potentially new centers after the swap.
            swapClusterAssignments.assignAll(data, centers);

            // The cost after the swap.
            double cost = swapClusterAssignments.calcCost();

            printf("Swaping cluster %3ld with point %3ld result in cost %0.5f\n", c, p, cost);

            if (cost < bestCost)
            {
                bestCost = cost;
                bestCenters = centers;
                std::cout << "Found new best centers: \n" << bestCenters << "\n";
            }
        }
    }
    
    return result;
}
