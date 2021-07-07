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

    auto rings = this->makeRings(clusterAssignments);

    addShortfallPointsToCoreset(clusterAssignments, rings, coresetPoints);

    auto groups = makeGroups(clusterAssignments, rings, 4);

    groupOvershotPoints(clusterAssignments, rings, groups);
    
    auto totalCost = clusterAssignments.getTotalCost();
    for (size_t i = 0; i < groups->size(); i++)
    {
        auto group = groups->get(i);
        auto groupPoints = group->getPoints();
        auto groupCost = group->calcTotalCost();
        auto normalizedGroupCost = groupCost / totalCost;
        auto numSamples = static_cast<size_t>(ceil(T * normalizedGroupCost));
        auto sampledPoints = random.choice(groupPoints, numSamples);

        printf("Group j=%ld l=%2d:   number of points=%2ld   cost=%2.4f   normalized cost=%0.4f   samples=%ld \n",
            group->RangeValue, group->RingRangeValue, groupPoints.size(), group->calcTotalCost(), normalizedGroupCost, 
            numSamples
        );

        printf("  Sampled points from group:\n");
        for (size_t i = 0; i < sampledPoints.size(); i++)
        {
            printf("    Point index %ld\n", sampledPoints[i]->PointIndex);
        }
    }

    // printPythonCodeForVisualisation(result, rings);
}

std::shared_ptr<RingSet>
GroupSampling::makeRings(const clustering::ClusterAssignmentList &clusterAssignments)
{
    const int ringRangeStart = -static_cast<int>(floor(std::log10(static_cast<double>(beta))));
    const int ringRangeEnd = -ringRangeStart;
    const auto n = clusterAssignments.getNumberOfPoints();
    const auto k = clusterAssignments.getNumberOfClusters();

    auto rings = std::make_shared<RingSet>(ringRangeStart, ringRangeEnd, k);

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
            auto ring = rings->findOrCreate(c, l, averageClusterCost);

            // Add point if cost(p, A) is within bounds i.e. between Δ_c*2^l and Δ_c*2^(l+1) 
            if (ring->tryAddPoint(p, costOfPoint))
            {
                pointPutInRing = true;

                // Since a point cannot belong to multiple rings, there is no need to
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
                // Track shortfall points: below l's lower range i.e. l<log⁡(1/β)
                rings->addShortfallPoint(p, c, costOfPoint, innerMostRingCost);
            }
            else if (costOfPoint > outerMostRingCost)
            {
                // Track overshot points: above l's upper range i.e., l>log⁡(β)
                rings->addOvershotPoint(p, c, costOfPoint, innerMostRingCost);
            }
            else
            {
                throw std::logic_error("Point should either belong to a ring or be ringless.");
            }
        }
    }

    return rings;
}

void GroupSampling::addShortfallPointsToCoreset(const clustering::ClusterAssignmentList &clusters, const std::shared_ptr<RingSet> rings, std::vector<WeightedPoint> &coresetPoints)
{
    // Handle points whose costs are below the lowest ring range i.e. l < log(1/beta).
    // These are called shortfall points because they fall short of being captured by the
    // inner-most ring. These points are snapped to the center of the assigned cluster by
    // adding the centers to the coreset weighted by the number of shortfall points of
    // that cluster.
    auto k = clusters.getNumberOfClusters();

    for (size_t c = 0; c < k; c++)
    {
        // The number of shortfall points for cluster `c`
        auto nShortfallPoints = rings->getNumberOfShortfallPoints(c);

        if (nShortfallPoints == 0)
        {
            // Not all clusters may have shortfall points so skip those.
            continue;
        }

        // The weight of the coreset point for the center of cluster `c`
        double weight = static_cast<double>(nShortfallPoints);

        // Add center to the coreset.
        printf("Add center of cluster %ld with weight %0.2f to coreset.\n", c, weight);
        coresetPoints.push_back(WeightedPoint(c, weight, true));
    }
}

void GroupSampling::groupOvershotPoints(const clustering::ClusterAssignmentList &clusters, const std::shared_ptr<RingSet> rings, std::shared_ptr<GroupSet> groups)
{
    size_t numberOfGroups = 5;
    auto k = clusters.getNumberOfClusters();
    double kDouble = static_cast<double>(k);
    double totalCost = rings->computeCostOfOvershotPoints();

    printf("\n\nGrouping overshot points, cost(O) = %0.5f\n", totalCost);
    
    for (size_t c = 0; c < k; c++)
    {
        double clusterCost = rings->computeCostOfOvershotPoints(c);
        auto points = rings->getOvershotPoints(c);
        
        printf("    Cluster i=%ld  - cost(C_i ⋂ O) = %0.4f     |C_i ⋂ O| = %ld\n", c, clusterCost, points.size());

        if (points.size() == 0)
        {
            // If no overshot points in the current cluster then go to next cluster.
            continue;
        }

        for (size_t j = 0; j < numberOfGroups; j++)
        {
            double jDouble = static_cast<double>(j);
            double lowerBound = 1/kDouble * pow(2, -jDouble    ) * totalCost;
            double upperBound = 1/kDouble * pow(2, -jDouble + 1) * totalCost;
            bool shouldAddPointsIntoGroup = false;

            if (j == 0)
            {
                // Group 0 has no upper bound. Notice this can be written as two-liners,
                // but is expanded to make it easier to read the code.
                shouldAddPointsIntoGroup = clusterCost >= lowerBound; 
                printf("\n      Group j=%ld    lowerBoundCost=%0.4f\n", j, lowerBound);
            }
            else
            {
                shouldAddPointsIntoGroup = clusterCost >= lowerBound && clusterCost < upperBound;
                printf("\n      Group j=%ld    lowerBoundCost=%0.4f   upperBoundCost=%0.4f\n", j, lowerBound, upperBound);
            }

            if (shouldAddPointsIntoGroup)
            {
                auto l = std::numeric_limits<int>().max();
                auto group = groups->create(j, l, lowerBound, upperBound);
                
                printf("            Adding %ld points to G[l=%d, j=%ld]\n", points.size(), l, j);

                for (size_t i = 0; i < points.size(); i++)
                {
                    auto point = points[i];
                    group->addPoint(point->PointIndex, point->ClusterIndex, point->PointCost);
                }
            }
        }
    }
}

std::shared_ptr<GroupSet>
GroupSampling::makeGroups(const clustering::ClusterAssignmentList &clusters, const std::shared_ptr<RingSet> rings, const size_t numberOfGroups)
{
    auto groups = std::make_shared<GroupSet>();
    auto k = static_cast<double>(clusters.getNumberOfClusters());
    for (int l = rings->RangeStart; l <= rings->RangeEnd; l++)
    {
        double ringCost = rings->calcRingCost(l);
        size_t nRingPointsForAllClusters = rings->countRingPoints(l);
        size_t nGroupedPoints = 0;
        printf("\n\nRing l=%d   -  cost(R_l) = %0.4f   -   |R_l| = %ld\n", l, ringCost, nRingPointsForAllClusters);
        
        for (size_t c = 0; c < k; c++)
        {
            auto ring = rings->find(c, l);
            auto clusterCost = ring->getTotalCost();
            auto ringPoints = ring->getPoints();

            printf("    Cluster i=%ld  - cost(R_{l,i}) = %0.4f     |R_{l,i}| = %ld\n", c, clusterCost, ring->countPoints());

            if (ring->countPoints() == 0)
            {
                // If nothing is captured by the current ring, then continue to the next ring.
                continue;
            }

            for (size_t j = 0; j < numberOfGroups; j++)
            {
                double jDouble = static_cast<double>(j);
                double lowerBound = 1/k * pow(2, -jDouble    ) * ringCost;
                double upperBound = 1/k * pow(2, -jDouble + 1) * ringCost;
                bool shouldAddPointsIntoGroup = false;

                if (j == 0)
                {
                    // Group 0 has no upper bound. Notice this can be written as two-liners,
                    // but is expanded to make it easier to read the code.
                    shouldAddPointsIntoGroup = clusterCost >= lowerBound; 
                    //printf("\n      Group j=%ld    lowerBoundCost=%0.4f\n", j, lowerBound, upperBound);
                }
                else
                {
                    shouldAddPointsIntoGroup = clusterCost >= lowerBound && clusterCost < upperBound;
                    //printf("\n      Group j=%ld    lowerBoundCost=%0.4f   upperBoundCost=%0.4f\n", j, lowerBound, upperBound);
                }
                
                if (shouldAddPointsIntoGroup)
                {
                    // Points which belong to cluster `c` and ring `l`
                    printf("            Adding %ld points to G[l=%d, j=%ld]\n", ringPoints.size(), l, j);

                    auto group = groups->create(j, l, lowerBound, upperBound);
                    for (size_t i = 0; i < ringPoints.size(); i++)
                    {
                        auto ringPoint = ringPoints[i];
                        group->addPoint(ringPoint->PointIndex, ringPoint->ClusterIndex, ringPoint->Cost);
                        nGroupedPoints++;
                    }
                }
            }
        }

        if (nRingPointsForAllClusters != nGroupedPoints)
        {
            printf("Not all points in ring l=%d are put in a group. ", l);
            printf("Number of points in the ring is %ld but only %ld points are grouped.\n",
                nRingPointsForAllClusters, nGroupedPoints);
        }
        assert(nRingPointsForAllClusters == nGroupedPoints);

    }

    return groups;
}

void
GroupSampling::printPythonCodeForVisualisation(std::shared_ptr<clustering::ClusteringResult> result, std::shared_ptr<RingSet> rings)
{
    auto clusterAssignments = result->getClusterAssignments();
    auto centers = result->getCentroids();
    auto k = clusterAssignments.getNumberOfClusters();
    auto n = clusterAssignments.getNumberOfPoints();

    printf("k = %ld\n", k);
    printf("cluster_labels = [");
    for (size_t p = 0; p < n; p++)
    {
        printf("%ld, ", clusterAssignments.getCluster(p));
    }
    printf("]\n");

    printf("cluster_centers = np.array([\n");
    for (size_t c = 0; c < centers.rows(); c++)
    {
        printf("  [");
        for (size_t d = 0; d < centers.columns(); d++)
        {
            printf("%0.5f, ", centers.at(c, d));
        }
        printf("],\n");
    }
    printf("])\n");
    
    printf("ring_ranges = [");
    for (int l = rings->RangeStart; l <= rings->RangeEnd; l++)
    {
        printf("%d, ", l);
    }
    printf("]\n");

    printf("rings = np.array([\n");
    for (size_t c = 0; c < k; c++)
    {
        printf("  [");
        for (int l = rings->RangeStart; l <= rings->RangeEnd; l++)
        {
            auto ring = rings->find(c, l);
            printf("%0.5f, ", ring->getLowerBoundCost());

        }
        printf("],\n");
    }
    printf("])\n");

    
    printf("rings_upper_bounds = np.array([\n");
    for (size_t c = 0; c < k; c++)
    {
        printf("  [");
        for (int l = rings->RangeStart; l <= rings->RangeEnd; l++)
        {
            auto ring = rings->find(c, l);
            printf("%0.5f, ", ring->getUpperBoundCost());

        }
        printf("],\n");
    }
    printf("])\n");
}
