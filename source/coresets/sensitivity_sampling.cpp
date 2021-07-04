#include <coresets/sensitivity_sampling.hpp>

using namespace coresets;

std::shared_ptr<blaze::DynamicVector<double>>
calcCenterWeights(
    const clustering::ClusterAssignmentList &clusterAssignments, const size_t targetCoresetPoints,
    std::shared_ptr<blaze::DynamicVector<size_t>> sampledIndices)
{
    auto sumOfCosts = clusterAssignments.getTotalCost();
    auto numberOfClusters = clusterAssignments.getNumberOfClusters();

    // Initialise an array to store center weights w_i
    auto centerWeights = std::make_shared<blaze::DynamicVector<double>>(numberOfClusters);
    centerWeights->reset();
    
    // For each of the T sampled points...
    for (auto &&p : *sampledIndices)
    {
        // Find point p's assigned cluster C_i.
        size_t clusterOfPointP = clusterAssignments.getCluster(p);

        // Find cost(p, A)
        double costPOfA = clusterAssignments.getPointCost(p);

        // Compute cost(A)/(T*cost(p,A))
        double weightContributionOfP = sumOfCosts / (targetCoresetPoints * costPOfA);

        printf("Point %3ld contributes %.5f to cluster %ld  ", p, weightContributionOfP, clusterOfPointP);

        // Sum it up: sum_{p sampled and p in C_i}   cost(A)/(T*cost(p,A))
        (*centerWeights)[clusterOfPointP] += weightContributionOfP;

        printf("  =>  w_%ld = %.5f\n", clusterOfPointP, (*centerWeights)[clusterOfPointP]);
    }

    printf("\n\n");

    // For each of the k' centers, compute the center weights.
    for (size_t c = 0; c < numberOfClusters; c++)
    {
        // Find precomputed center weight: w_i
        double w_i = (*centerWeights)[c];

        // Compute |C_i|
        size_t numberOfPointsInCluster = clusterAssignments.countPointsInCluster(c);

        // Compute max(0, |C_i| - w_i)
        double centerWeight = blaze::max(0.0, static_cast<double>(numberOfPointsInCluster) - w_i);

        // Update the center weight.
        (*centerWeights)[c] = centerWeight;

        printf("|C_%ld| = %3ld,  w_%ld = %2.5f,  new w_%ld = %2.5f \n", c, numberOfPointsInCluster, c, w_i, c, centerWeight);
    }

    return centerWeights;
}

CoresetResult::CoresetResult()
{
}

SensitivitySampling::SensitivitySampling(const int randomSeed) : random(randomSeed)
{
}


std::vector<WeightedPoint>
SensitivitySampling::calcCoresetPoints(const std::shared_ptr<clustering::ClusteringResult> result, const size_t targetCoresetPoints)
{
    // Generated coresets points that the method returns.
    std::vector<WeightedPoint> coresetPoints;

    // Number of points that should be in the coreset.
    //auto T = targetCoresetPoints;
    auto centers = result->getCentroids();
    auto clusterAssignments = result->getClusterAssignments();
    auto n = clusterAssignments.getNumberOfPoints();

    // Step 2a: compute cost(p, A). The cost of each point is the
    // distance between the point and its closest centroid.
    auto costs = clusterAssignments.getCentroidDistances();

    // Step 2b: compute cost(A). Assume it is the sum of all costs.
    auto sumOfCosts = clusterAssignments.getTotalCost();

    // Step 2c: compute the sampling distribution: cost(p, A)/cost(A)
    auto samplingDistribution = costs / sumOfCosts;

    // TODO: Investigate why small weights generate samples that are all zeros.
    auto sampledIndices = random.choice(targetCoresetPoints, n, samplingDistribution * 100);
    std::cout << "Sampled T points: \n"
              << (*sampledIndices) << "\n";

    // Loop through the sampled points and calculate
    // the weight associated with each of these points.
    for (size_t j = 0; j < targetCoresetPoints; j++)
    {
        size_t sampledPointIndex = (*sampledIndices)[j];

        // We scale the cost of the sampled point by a factor of T i.e. T * cost(p,A)
        double scaledCostPofA = targetCoresetPoints * costs[sampledPointIndex];

        // The weight of the sampled point is now: cost(A) / (T*cost(p,A))
        double weight = sumOfCosts / scaledCostPofA;

        coresetPoints.push_back(WeightedPoint(sampledPointIndex, weight, false));

        printf("Sampled point %3ld gets weight %.5f \n", sampledPointIndex, weight);
    }

    auto numberOfClusters = clusterAssignments.getNumberOfClusters();
    auto centerWeights = calcCenterWeights(clusterAssignments, targetCoresetPoints, sampledIndices);

    for (size_t c = 0; c < numberOfClusters; c++)
    {
        auto weight = (*centerWeights)[c];
        coresetPoints.push_back(WeightedPoint(c, weight, true));
    }

    return coresetPoints;
}

std::shared_ptr<CoresetResult>
SensitivitySampling::run(const blaze::DynamicMatrix<double> &data)
{
    uint k = 3;
    uint kPrime = 2 * k;
    uint targetCoresetPoints = 20; // T is the number of sampled points. It is hyperparam. Usually T=200*k

    // Step 1: Run k-means++ to get the initial solution A.
    clustering::KMeans kMeansAlg(kPrime, true, 100U, 0.0001, 42);
    auto result = kMeansAlg.run(data);

    auto coresetPoints = calcCoresetPoints(result, targetCoresetPoints);

    return std::make_shared<CoresetResult>();
}
