#include <coresets/group_sampling.hpp>

using namespace coresets;

GroupSampling::GroupSampling() : beta(200)
{
}

void GroupSampling::run(const blaze::DynamicMatrix<double> &data)
{
    utils::Random random;
    const uint kPrime = 3;
    const uint k = kPrime; // TODO: Should be k = 2 * kPrime;
    const uint T = 20;     // T is the number of sampled points. It is hyperparam. Usually T=200*k
    const size_t n = data.rows();
    const size_t d = data.columns();
    auto coreset = std::make_shared<Coreset>(T, d);
    std::vector<WeightedPoint> coresetPoints;

    // Step 1: Run k-means++ to get the initial solution A.
    clustering::KMeans kMeansAlg(k);
    auto result = kMeansAlg.run(data);
    auto centers = result->getCentroids();

    auto clusterAssignments = result->getClusterAssignments();

    // Step 2: Compute the average cost for each cluster.
    auto averageClusterCosts = clusterAssignments.calcAverageClusterCosts();

    for (size_t c = 0; c < k; c++)
    {
        printf("Cluster %ld's average cost is %.4f\n", c, (*averageClusterCosts)[c]);
    }

    auto rings = this->makeRings(result);

    addInnerMostRingPoints(clusterAssignments, rings, coresetPoints);

    addOuterMostRingPoints(data, rings, coresetPoints);

    auto groups = makeGroups(clusterAssignments, rings, 4);
    auto totalCost = clusterAssignments.getTotalCost();
    for (size_t i = 0; i < groups->size(); i++)
    {
        auto group = groups->get(i);
        auto groupPoints = group->getPoints();
        auto groupCost = group->calcTotalCost();
        auto normalizedGroupCost = groupCost / totalCost;
        auto numSamples = static_cast<size_t>(ceil(T * normalizedGroupCost));

        printf("Group j=%ld l=%2d:   number of points=%2ld   cost=%2.4f   normalized cost=%0.4f   samples=%ld \n",
            group->RangeValue, group->RingRangeValue, groupPoints.size(), group->calcTotalCost(), normalizedGroupCost, 
            numSamples
        );
    }
}

std::shared_ptr<RingSet>
GroupSampling::makeRings(const std::shared_ptr<clustering::ClusteringResult> clusters)
{
    auto clusterAssignments = clusters->getClusterAssignments();
    const int ringRangeStart = -static_cast<int>(floor(std::log10(static_cast<double>(beta))));
    const int ringRangeEnd = -ringRangeStart;
    auto rings = std::make_shared<RingSet>(ringRangeStart, ringRangeEnd);
    const auto n = clusterAssignments.getNumberOfPoints();
    const auto k = clusterAssignments.getNumberOfClusters();

    // Step 2: Compute the average cost for each cluster.
    auto averageClusterCosts = clusterAssignments.calcAverageClusterCosts();

    for (size_t p = 0; p < n; p++)
    {
        // The cluster index of point `p`
        size_t c = clusterAssignments.getCluster(p);

        // The cost of point `p`: cost(p, A)
        double costOfPoint = clusterAssignments.getPointCost(p);

        // The average cost of cluster `c`: Δ_c
        double averageClusterCost = (*averageClusterCosts)[c];

        bool pointPutInRing = false;
        for (int l = ringRangeStart; l <= ringRangeEnd; l++)
        {
            // Ring upper bound := Δ_c * 2^l
            double ringLowerBound = averageClusterCost * std::pow(2, l);

            // Ring upper bound := Δ_c * 2^(l+1)
            double ringUpperBound = averageClusterCost * std::pow(2, l + 1);

            // If cost(p, A) is between Δ_c*2^l and Δ_c*2^(l+1) ...
            if (costOfPoint >= ringLowerBound && costOfPoint < ringUpperBound)
            {
                pointPutInRing = true;
                auto ring = std::make_shared<InternalRing>(p, c, costOfPoint, l, ringLowerBound, ringUpperBound);
                rings->add(ring);

                // Since a point cannot belong to multiple rings, there is no need to look
                // test whether the point `p` falls within the ring of the next range l+1.
                break;
            }
        }

        if (pointPutInRing == false)
        {
            double innerMostRingCost = averageClusterCost * std::pow(2, ringRangeStart);
            double outerMostRingCost = averageClusterCost * std::pow(2, ringRangeEnd + 1);

            if (costOfPoint < innerMostRingCost)
            {
                // Step 5: Handle points below l's lower range i.e. l<log⁡(1/β)
                auto ring = std::make_shared<ExternalRing>(p, c, costOfPoint, innerMostRingCost, true);
                rings->add(ring);
            }
            else if (costOfPoint > outerMostRingCost)
            {
                // Step 6: Handle points above l's upper range i.e., l>log⁡(β)
                auto ring = std::make_shared<ExternalRing>(p, c, costOfPoint, outerMostRingCost, false);
                rings->add(ring);
            }
            else
            {
                throw std::logic_error("Point does not belong to internal nor external ring. Program logic error.");
            }
        }
    }

    return rings;
}

void GroupSampling::addInnerMostRingPoints(const clustering::ClusterAssignmentList &clusters, std::shared_ptr<RingSet> rings, std::vector<WeightedPoint> &coresetPoints)
{
    // Step 5: Handle points that are below the lowest ring range i.e. l < log(1/beta).
    // These are called inner-most external rings. Since these points are very close to
    // their corresponding cluster centers, we add their centers to the coreset weighted
    // by the number of points in the inner-most external ring of that cluster.
    auto k = clusters.getNumberOfClusters();

    for (size_t c = 0; c < k; c++)
    {
        // The number of points in the inner-most external ring for cluster `c`
        auto innerRingPointClusterCount = rings->getNumberOfInnerRingPoints(c);

        if (innerRingPointClusterCount == 0)
        {
            // Not all clusters may have points in the inner-most external ring, so
            // do not add center of the curresponding cluster `c` to the coreset.
            continue;
        }

        // The weight of the coreset point for the center of cluster `c`
        double weight = static_cast<double>(innerRingPointClusterCount);

        // Add center to the coreset.
        printf("Add center of cluster %ld with weight %0.2f to coreset.\n", c, weight);
        coresetPoints.push_back(WeightedPoint(c, weight, true));
    }
}

void GroupSampling::addOuterMostRingPoints(const blaze::DynamicMatrix<double> &data, std::shared_ptr<RingSet> rings, std::vector<WeightedPoint> &coresetPoints)
{
    // First, make a matrix that consists of points in the dataset that are captured by the outer-most ring.
    std::vector<size_t> pointIndicesOutsideAllRings = rings->getPointsOutsideAllRings();
    auto pPrimeSize = pointIndicesOutsideAllRings.size();
    if (pPrimeSize > 0)
    {
        blaze::DynamicMatrix<double> pointsOutsideAllRings(pPrimeSize, data.columns());
        std::cout << "Points outside all rings = [ ";
        for (size_t i = 0; i < pPrimeSize; i++)
        {
            size_t p = pointIndicesOutsideAllRings[i];
            blaze::row(pointsOutsideAllRings, i) = blaze::row(data, p);
            std::cout << p << ", ";
        }
        std::cout << "]\n";
        
        std::cout << "Points: \n" << pointsOutsideAllRings << "\n";

        // TODO: Run sensitivity sampling.
        // TODO: Add results of Sensitivity Sampling into coresetPoints
    }
}

std::shared_ptr<GroupSet>
GroupSampling::makeGroups(const clustering::ClusterAssignmentList &clusters, const std::shared_ptr<RingSet> rings, const size_t numberOfGroups)
{
    auto groups = std::make_shared<GroupSet>();
    auto k = static_cast<double>(clusters.getNumberOfClusters());
    auto clusterCosts = clusters.calcAverageClusterCosts();
    for (int l = rings->RangeStart; l <= rings->RangeEnd; l++)
    {
        double ringCost = rings->calcRingCost(l);
        printf("\n\nRing cost for l=%d is %0.4f\n", l, ringCost);
        
        for (size_t j = 0; j < numberOfGroups; j++)
        {
            double jDouble = static_cast<double>(j);
            double lowerBound = 1/k * pow(2, -jDouble    ) * ringCost;
            double upperBound = 1/k * pow(2, -jDouble + 1) * ringCost;

            printf("   Group j=%ld -> lowerBoundCost=%0.4f  upperBoundCost=%0.4f\n", j, lowerBound, upperBound);

            for (size_t c = 0; c < k; c++)
            {
                auto clusterCost = (*clusterCosts)[c];
                
                printf("        Cluster C_%ld: cost(C_%ld)=%0.4f\n", c, c, clusterCost);
                if (clusterCost >= lowerBound && clusterCost < upperBound)
                {
                    // Points which belong to cluster `c` and ring `l`
                    auto ringPoints = rings->getInternalRingPointsInCluster(c, l);

                    printf("            Cluster C_%ld has %ld ring points for range %d\n", c, ringPoints.size(), l);

                    if (ringPoints.size() > 0)
                    {
                        auto group = groups->create(j, l, lowerBound, upperBound);
                        for (size_t i = 0; i < ringPoints.size(); i++)
                        {
                            auto ringPoint = ringPoints[i];
                            group->addPoint(ringPoint->PointIndex, ringPoint->ClusterIndex, ringPoint->PointCost);
                        }
                    }
                }
            }
        }
    }

    return groups;
}
