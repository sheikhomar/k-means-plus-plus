#include <clustering/kmeans.hpp>

using namespace clustering;

KMeans::KMeans(uint k, bool kpp, uint miter, double convDiff, int randSeed) : numOfClusters(k), initKMeansPlusPlus(kpp), maxIterations(miter), convergenceDiff(convDiff), random(randSeed)
{
}

std::shared_ptr<ClusteringResult>
KMeans::run(const blaze::DynamicMatrix<double> &data)
{
  auto centroids = (this->initKMeansPlusPlus ? this->initCentroidsKMeansPlusPlus(data) : this->initCentroidsNaive(data));
  return this->runLloydsAlgorithm(data, centroids);
}

blaze::DynamicMatrix<double>
KMeans::initCentroidsNaive(const blaze::DynamicMatrix<double> &dataPoints)
{
  auto n = dataPoints.rows();
  auto d = dataPoints.columns();
  auto k = this->numOfClusters;
  auto randomPointGenerator = random.getIndexer(n);

  blaze::DynamicMatrix<double> centers(k, d);

  for (size_t c = 0; c < k; c++)
  {
    // Pick a random point p as the cluster center.
    auto randomPoint = randomPointGenerator.next();
    blaze::row(centers, c) = blaze::row(dataPoints, randomPoint);
  }

  return centers;
}

blaze::DynamicMatrix<double>
KMeans::initCentroidsKMeansPlusPlus(const blaze::DynamicMatrix<double> &matrix)
{
  size_t n = matrix.rows();
  size_t d = matrix.columns();
  auto k = this->numOfClusters;
  auto randomPointGenerator = random.getIndexer(n);

  blaze::DynamicMatrix<double> centroids(k, d);
  std::vector<uint> availableIndices(n);

  // Fill with 0, 1, 2, ..., N.
  std::iota(availableIndices.begin(), availableIndices.end(), 0);

  for (size_t c = 0; c < k; c++)
  {
    size_t centroidIndex = 0;

    if (c == 0)
    {
      // Pick the first centroid uniformly at random.
      centroidIndex = randomPointGenerator.next();
    }
    else
    {

      blaze::DynamicVector<double> smallestDistances(n);

      // For each point, find the distance to the nearest centroid for all
      // the centroids that are select so far.
      for (uint p : availableIndices)
      {
        double smallestDistance = std::numeric_limits<double>::max();

        // Loop through previously selected clusters.
        for (size_t c2 = 0; c2 < c; c2++)
        {

          // Compute the L2 norm between point p and centroid c2.
          const double distance = blaze::norm(blaze::row(matrix, p) - blaze::row(centroids, c2));

          // Decide if current distance is better.
          if (distance < smallestDistance)
          {
            smallestDistance = distance;
          }
        }

        smallestDistances[p] = smallestDistance;
      }

      // Pick a point based on a weighted probability

      // Square distances
      smallestDistances *= smallestDistances;

      // Normalise.
      smallestDistances /= blaze::sum(smallestDistances);

      // Pick the index of a point randomly selected based on the weights.
      size_t nextClusterCandidate = random.choice(smallestDistances);

      // Assign centroid index.
      centroidIndex = availableIndices[nextClusterCandidate];
    }

    std::cout << "Centroid index for " << c << " => " << centroidIndex << "\n";

    // Copy point over centroids matrix.
    blaze::row(centroids, c) = blaze::row(matrix, centroidIndex);

    // Remove it from the candidate list so the point cannot be
    // picked as another centroid.
    boost::remove_erase(availableIndices, centroidIndex);
  }

  return centroids;
}

std::shared_ptr<ClusteringResult>
KMeans::runLloydsAlgorithm(const blaze::DynamicMatrix<double> &matrix, blaze::DynamicMatrix<double> centroids)
{
  size_t n = matrix.rows();
  auto k = this->numOfClusters;

  blaze::DynamicVector<size_t> clusterMemberCounts(k);
  ClusterAssignmentList cal(n, this->numOfClusters);
  
  for (size_t i = 0; i < this->maxIterations; i++)
  {
    // For each data point, assign the centroid that is closest to it.
    for (size_t p = 0; p < n; p++)
    {
      double bestDistance = std::numeric_limits<double>::max();
      size_t bestCluster = 0;

      // Loop through all the clusters.
      for (size_t c = 0; c < k; c++)
      {

        // Compute the L2 norm between point p and centroid c.
        const double distance = blaze::norm(blaze::row(matrix, p) - blaze::row(centroids, c));

        // Decide if current distance is better.
        if (distance < bestDistance)
        {
          bestDistance = distance;
          bestCluster = c;
        }
      }

      // Assign cluster to the point p.
      cal.assign(p, bestCluster, bestDistance);
    }

    // Move centroids based on the cluster assignments.

    // First, save a copy of the centroids matrix.
    blaze::DynamicMatrix<double> oldCentrioids(centroids);

    // Set all elements to zero.
    centroids = 0;           // Reset centroids.
    clusterMemberCounts = 0; // Reset cluster member counts.

    for (size_t p = 0; p < n; p++)
    {
      const size_t c = cal.getCluster(p);
      blaze::row(centroids, c) += blaze::row(matrix, p);
      clusterMemberCounts[c] += 1;
    }

    for (size_t c = 0; c < k; c++)
    {
      const auto count = std::max<size_t>(1, clusterMemberCounts[c]);
      blaze::row(centroids, c) /= count;
    }

    std::cout << "Centroids after iteration " << i << ": \n"
              << centroids << "\n";

    // Compute the Frobenius norm
    auto diffAbsMatrix = blaze::abs(centroids - oldCentrioids);
    auto diffAbsSquaredMatrix = blaze::pow(diffAbsMatrix, 2); // Square each element.
    auto frobeniusNormDiff = blaze::sqrt(blaze::sum(diffAbsSquaredMatrix));

    std::cout << "Frobenius norm of centroids difference: " << frobeniusNormDiff << "!\n";

    if (frobeniusNormDiff < this->convergenceDiff)
    {
      std::cout << "Stopping k-Means as centroids do not improve. Frobenius norm Diff: " << frobeniusNormDiff << "\n";
      break;
    }
  }

  return std::make_shared<ClusteringResult>(cal, centroids);
}
