#include <kmeans/clustering_result.hpp>

using namespace kmeans;

ClusteringResult::ClusteringResult(const kmeans::ClusterAssignmentList &assignments, blaze::DynamicMatrix<double> &finalCentroids) :
    clusterAssignments(assignments), centroids(finalCentroids)
{
}

const kmeans::ClusterAssignmentList&
ClusteringResult::getClusterAssignments()
{
    return this->clusterAssignments;
}

const blaze::DynamicMatrix<double>&
ClusteringResult::getCentroids()
{
    return this->centroids;
}
