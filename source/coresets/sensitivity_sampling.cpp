#include <coresets/sensitivity_sampling.hpp>

using namespace coresets;

CoresetResult::CoresetResult()
{
}

SensitivySampling::SensitivySampling()
{
}

std::shared_ptr<CoresetResult>
SensitivySampling::run(const blaze::DynamicMatrix<double> &data)
{
    utils::Random random(42);
    uint k = 3;
    uint kPrime = 2 * k;
    uint T = 20;         // T is the number of sampled points. It is hyperparam. Usually T=200*k
    size_t n = data.rows();

    // Step 1: Run k-means++ to get the initial solution A.
    kmeans::KMeans kMeansAlg(kPrime, true, 100U, 0.0001, 42);
    auto result = kMeansAlg.run(data);
    auto A = result->getCentroids();

    auto clusterAssignments = result->getClusterAssignments();
    auto numOfPoints = clusterAssignments.getNumberOfPoints();

    // Step 2a: compute cost(p, A). The cost of each point is the
    // distance between the point and its closest centroid.
    auto costs = clusterAssignments.getCentroidDistances();

    // Step 2b: compute cost(A). Assume it is the sum of all costs.
    auto sumOfCosts = blaze::sum(costs);

    // Step 2c: compute the sampling distribution: cost(p, A)/cost(A)
    auto samplingDistribution = costs / sumOfCosts;

    // TODO: Investigate why small weights generate samples that are all zeros.
    auto sampledIndices = random.choice(T, n, samplingDistribution * 100);
    std::cout << "Sampled T points: \n"
              << (*sampledIndices) << "\n";

    // Initialise an array to store the weights of the sampled points.
    blaze::DynamicVector<double> sampledPointWeights(T);

    for (size_t j = 0; j < T; j++)
    {
        size_t sampledPointIndex = (*sampledIndices)[j];

        // We scale the cost of the sampled point by a factor of T i.e. T * cost(p,A)
        double scaledCostPofA = T * costs[sampledPointIndex];

        // The weight of the sampled point is now: cost(A) / (T*cost(p,A))
        sampledPointWeights[j] = sumOfCosts / scaledCostPofA;
    }

    // Initialise an array to store center weights w_i
    blaze::DynamicVector<double> centerWeights(kPrime);

    // For each of the T sampled points...
    for (auto &&p : *sampledIndices)
    {
        // Find point p's assigned cluster C_i.
        size_t clusterOfPointP = clusterAssignments.getCluster(p);

        // Find cost(p, A)
        double costPOfA = costs[p];

        // Compute cost(A)/(T*cost(p,A))
        double weightContributionOfP = sumOfCosts / (T * costPOfA);

        std::cout << "Point " << p << " contributes " << weightContributionOfP << " to cluster " << clusterOfPointP << "\n";

        // Sum it up: sum_{p sampled and p in C_i}   cost(A)/(T*cost(p,A))
        centerWeights[clusterOfPointP] += weightContributionOfP;

        std::cout << "  w_" << clusterOfPointP << " = " << centerWeights[clusterOfPointP] << "\n";
    }

    // For each of the k' centers, compute the center weights.
    for (size_t c = 0; c < kPrime; c++)
    {
        // Find precomputed center weight: w_i
        double w_i = centerWeights[c];

        // Compute |C_i|
        size_t numberOfPointsInCluster = clusterAssignments.countPointsInCluster(c);

        // Compute max(0, |C_i| - w_i)
        double centerWeight = blaze::max(0, numberOfPointsInCluster - w_i);

        // Update the center weight.
        centerWeights[c] = centerWeight;
    }

    return std::make_shared<CoresetResult>();
}
