#include <kmeans/cluster_assignment_list.hpp>

using namespace kmeans;

ClusterAssignmentList::ClusterAssignmentList(uint n, uint k) : numOfPoints(n), numOfClusters(k), clusters(n), distances(n)
{
}

void ClusterAssignmentList::assign(size_t pointIndex, size_t clusterIndex, double distance)
{
    // TODO: Ensure arguments are not out of range to avoid runtime errors.
    clusters[pointIndex] = clusterIndex;
    distances[pointIndex] = distance;
}

size_t
ClusterAssignmentList::getCluster(size_t pointIndex)
{
    // TODO: Ensure arguments are not out of range to avoid runtime errors.
    return clusters[pointIndex];
}

uint
ClusterAssignmentList::getNumberOfPoints()
{
    return this->numOfPoints;
}

uint
ClusterAssignmentList::getNumberOfClusters()
{
    return this->numOfClusters;
}
