#include <clustering/cluster_assignment_list.hpp>

using namespace clustering;

ClusterAssignmentList::ClusterAssignmentList(uint n, uint k) : numOfPoints(n), numOfClusters(k), clusters(n), distances(n)
{
}

ClusterAssignmentList::ClusterAssignmentList(const ClusterAssignmentList& other) : 
    numOfPoints(other.numOfPoints), numOfClusters(other.numOfClusters), clusters(other.clusters), distances(other.distances)
{
}

void ClusterAssignmentList::assign(size_t pointIndex, size_t clusterIndex, double distance)
{
    // TODO: Ensure arguments are not out of range to avoid runtime errors.
    clusters[pointIndex] = clusterIndex;
    distances[pointIndex] = distance;
}

void
ClusterAssignmentList::assignAll(const blaze::DynamicMatrix<double> &dataPoints, const blaze::DynamicMatrix<double> &centers)
{
    auto n = this->numOfPoints;
    auto k = this->numOfClusters;

    // For each data point, assign the centroid that is closest to it.
    for (size_t p = 0; p < n; p++)
    {
      double bestDistance = std::numeric_limits<double>::max();
      size_t bestCluster = 0;

      // Loop through all the clusters.
      for (size_t c = 0; c < k; c++)
      {
        // Compute the L2 norm between point p and centroid c.
        const double distance = blaze::norm(blaze::row(dataPoints, p) - blaze::row(centers, c));

        // Decide if current distance is better.
        if (distance < bestDistance)
        {
          bestDistance = distance;
          bestCluster = c;
        }
      }

      // Assign cluster to the point p.
      this->assign(p, bestCluster, bestDistance);
    }
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

blaze::DynamicVector<double>&
ClusterAssignmentList::getCentroidDistances()
{
    return this->distances;
}

size_t
ClusterAssignmentList::countPointsInCluster(size_t clusterIndex)
{
    size_t count = 0;
    for (size_t p = 0; p < this->numOfPoints; p++)
    {
        if (clusters[p] == clusterIndex)
        {
            count++;
        }
    }
    
    return count;
}

double
ClusterAssignmentList::calcCost()
{
    return blaze::sum(this->distances);
}
