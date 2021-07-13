#include <clustering/kmeans.hpp>

using namespace clustering;

KMeans::KMeans(uint k, bool kpp, uint miter, double convDiff) : numOfClusters(k), initKMeansPlusPlus(kpp), maxIterations(miter), convergenceDiff(convDiff)
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
  utils::Random random;
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
KMeans::initCentroidsKMeansPlusPlus2(const blaze::DynamicMatrix<double> &matrix)
{
  utils::Random random;
  size_t n = matrix.rows();
  size_t d = matrix.columns();
  auto k = this->numOfClusters;

  // Compute squared pairwise distances between points.
  auto M = matrix * blaze::trans(matrix);
  blaze::DynamicVector<double> diagM(n);
  diagM = blaze::diagonal(M);
  blaze::DynamicVector<double> ones(n);
  ones = 1;
  auto h = diagM * blaze::trans(ones);
  blaze::DynamicMatrix<double> squaredPairwiseDistances(h + blaze::trans(h) - 2 * M);

  // Track which points are picked as centers.
  std::vector<size_t> pickedPointsAsCenters;
  pickedPointsAsCenters.reserve(k);

  for (size_t c = 0; c < k; c++)
  {
    size_t centerIndex = 0;

    if (c == 0)
    {
      // Pick the first centroid uniformly at random.
      auto randomPointGenerator = random.getIndexer(n);
      centerIndex = randomPointGenerator.next();
    }
    else
    {
      blaze::DynamicVector<double> weights(n);
      for (size_t p1 = 0; p1 < n; p1++)
      {
        double smallestDistance = std::numeric_limits<double>::max();

        // Loop through previously selected center points.
        for (size_t p2 : pickedPointsAsCenters)
        {
          // Notice that for points previously picked as centers, their
          // distances will be zero because the diagonal elements of the
          // pairwise distance matrix are all zeros: D[i,i] = 0.
          double distance = squaredPairwiseDistances.at(p1, p2);

          // Decide if current distance is better.
          if (distance < smallestDistance)
          {
            smallestDistance = distance;
          }
        }

        // Set the weight of a given point to be the smallest distance
        // to any of the previously selected center points. We want to
        // select points randomly such that points that are far from
        // any of the selected center points have higher likelihood of
        // being picked as the next candidate center.
        weights[p1] = smallestDistance;
      }

      // Normalise the weights.
      weights /= blaze::sum(weights);

      // Pick the index of a point randomly selected based on the weights.
      centerIndex = random.choice(weights);
    }

    std::cout << "Center index for " << c << " => " << centerIndex << "\n";
    pickedPointsAsCenters.push_back(centerIndex);
  }

  blaze::DynamicMatrix<double> centroids(k, d);
  for (size_t c = 0; c < k; c++)
  {
    size_t pointIndex = pickedPointsAsCenters[c];
    blaze::row(centroids, c) = blaze::row(matrix, pointIndex);
  }
  return centroids;
}

blaze::DynamicMatrix<double>
KMeans::initCentroidsKMeansPlusPlus(const blaze::DynamicMatrix<double> &matrix)
{
  utils::Random random;
  size_t n = matrix.rows();
  size_t d = matrix.columns();
  auto k = this->numOfClusters;

  // Track which points are picked as centers.
  std::vector<size_t> pickedPointsAsCenters;
  pickedPointsAsCenters.reserve(k);

  for (size_t c = 0; c < k; c++)
  {
    size_t centerIndex = 0;

    if (c == 0)
    {
      // Pick the first centroid uniformly at random.
      auto randomPointGenerator = random.getIndexer(n);
      centerIndex = randomPointGenerator.next();
    }
    else
    {
      blaze::DynamicVector<double> weights(n);
      for (size_t p1 = 0; p1 < n; p1++)
      {
        double smallestDistance = std::numeric_limits<double>::max();

        // Loop through previously selected clusters.
        for (size_t p2 : pickedPointsAsCenters)
        {
          double distance = 0.0;

          if (p1 != p2)
          {
            // Compute the squared L2 norm between point p1 and center p2.
            distance = blaze::sqrNorm(blaze::row(matrix, p1) - blaze::row(matrix, p2));
          }

          // Decide if current distance is better.
          if (distance < smallestDistance)
          {
            smallestDistance = distance;
          }
        }

        weights[p1] = smallestDistance;
      }

      // Normalise the weights.
      weights /= blaze::sum(weights);

      // Pick the index of a point randomly selected based on the weights.
      centerIndex = random.choice(weights);
    }

    std::cout << "Centroid index for " << c << " => " << centerIndex << "\n";
    pickedPointsAsCenters.push_back(centerIndex);
  }

  blaze::DynamicMatrix<double> centroids(k, d);
  for (size_t c = 0; c < k; c++)
  {
    size_t pointIndex = pickedPointsAsCenters[c];
    blaze::row(centroids, c) = blaze::row(matrix, pointIndex);
  }
  return centroids;
}

std::vector<size_t>
KMeans::pickInitialCentersViaKMeansPlusPlus(const blaze::DynamicMatrix<double> &matrix, const bool usePrecomputeDistances)
{
  utils::Random random;
  size_t n = matrix.rows();
  size_t d = matrix.columns();
  auto k = this->numOfClusters;

  blaze::DynamicMatrix<double> pairwiseDist;

  if (usePrecomputeDistances)
  {
    // Compute squared pairwise distances between points.
    printf("Precomputing pairwise squared distances between all points!\n");

    auto M = matrix * blaze::trans(matrix);
    blaze::DynamicVector<double> diagM(n);
    diagM = blaze::diagonal(M);
    blaze::DynamicVector<double> ones(n);
    ones = 1;
    auto h = diagM * blaze::trans(ones);
    pairwiseDist = h + blaze::trans(h) - 2 * M;
  }

  // Lambda function computes the squared L2 distance between any pair of points.
  // The function will automatically use any precomputed distance if it exists.
  auto calcSquaredL2Norm = [matrix, pairwiseDist, usePrecomputeDistances](size_t p1, size_t p2) -> double
  {
    if (p1 == p2)
    {
      return 0.0;
    }

    if (usePrecomputeDistances)
    {
      return pairwiseDist.at(p1, p2);
    }

    return blaze::sqrNorm(blaze::row(matrix, p1) - blaze::row(matrix, p2));
  };

  // Track which points are picked as centers.
  std::vector<size_t> pickedPointsAsCenters;
  pickedPointsAsCenters.reserve(k);

  for (size_t c = 0; c < k; c++)
  {
    size_t centerIndex = 0;

    if (c == 0)
    {
      // Pick the first centroid uniformly at random.
      auto randomPointGenerator = random.getIndexer(n);
      centerIndex = randomPointGenerator.next();
    }
    else
    {
      blaze::DynamicVector<double> weights(n);
      for (size_t p1 = 0; p1 < n; p1++)
      {
        double smallestDistance = std::numeric_limits<double>::max();

        // Loop through previously selected clusters.
        for (size_t p2 : pickedPointsAsCenters)
        {
          // Notice that for points previously picked as centers, their
          // distances will be zero because the diagonal elements of the
          // pairwise distance matrix are all zeros: D^2[i,i] = 0.
          double distance = calcSquaredL2Norm(p1, p2);

          // Decide if current distance is better.
          if (distance < smallestDistance)
          {
            smallestDistance = distance;
          }
        }

        // Set the weight of a given point to be the smallest distance
        // to any of the previously selected center points. We want to
        // select points randomly such that points that are far from
        // any of the selected center points have higher likelihood of
        // being picked as the next candidate center.
        weights[p1] = smallestDistance;
      }

      // Normalise the weights.
      weights /= blaze::sum(weights);

      // Pick the index of a point randomly selected based on the weights.
      centerIndex = random.choice(weights);
    }

    std::cout << "Center index for " << c << " => " << centerIndex << "\n";
    pickedPointsAsCenters.push_back(centerIndex);
  }

  return pickedPointsAsCenters;
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
