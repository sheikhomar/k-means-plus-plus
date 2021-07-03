#include <coresets/group_sampling.hpp>

using namespace coresets;

GroupSampling::GroupSampling()
{
}

void
GroupSampling::run(const blaze::DynamicMatrix<double> &data)
{
    utils::Random random(42);
    uint kPrime = 3;
    uint k = kPrime; // TODO: Should be k = 2 * kPrime;
    uint T = 20;         // T is the number of sampled points. It is hyperparam. Usually T=200*k
    size_t n = data.rows();
    size_t beta = 100;
    int ringRangeStart = -std::log10(static_cast<double>(beta));
    int ringRangeEnd = -ringRangeStart;

    // Step 1: Run k-means++ to get the initial solution A.
    clustering::KMeans kMeansAlg(k, true, 100U, 0.0001, 42);
    auto result = kMeansAlg.run(data);
    auto A = result->getCentroids();

    auto clusterAssignments = result->getClusterAssignments();

    // Step 2: Compute the average cost for each cluster.
    auto averageClusterCosts = clusterAssignments.calcAverageClusterCosts();

    for (size_t c = 0; c < k; c++)
    {
        printf("Cluster %ld's average cost is %.4f\n", c, (*averageClusterCosts)[c]);
    }

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
            double ringUpperBound = averageClusterCost * std::pow(2, l+1);

            // If cost(p, A) is between Δ_c*2^l and Δ_c*2^(l+1) ...
            if (costOfPoint >= ringLowerBound && costOfPoint < ringUpperBound) 
            {
                printf("Point %3ld with cost(p, A) = %0.4f  ->  R[%2d, %ld]  \n", p, costOfPoint, l, c);
                pointPutInRing = true;
                break; // A point cannot belong to multiple rings.
            }
        }

        if (pointPutInRing == false)
        {
            //printf("Point %3ld (cost=%0.5f) did not get be put in a ring: ", p, costOfPoint);
            printf("Point %3ld with cost(p, A) = %0.4f  -> no ring because ", p, costOfPoint);

            double innerMostRingCost = averageClusterCost * std::pow(2, ringRangeStart);
            double outerMostRingCost = averageClusterCost * std::pow(2, ringRangeEnd+1);

            if (costOfPoint < innerMostRingCost)
            {
                // Step 5: Handle points below l's lower range i.e. l<log⁡(1/β)

                printf(" the cost(p, A) falls below the cost range of inner most ring (%.4f)", innerMostRingCost);
            }
            else if (costOfPoint > outerMostRingCost)
            {
                // Step 6: Handle points above l's upper range i.e., l>log⁡(β)

                printf(" the cost(p, A) is above the cost range of outer most ring (%.4f)", outerMostRingCost);
            } 
            else
            {
                printf(" something went wrong!"); // TODO: Raise exception.
            }

            printf("\n");
        }
    }
}
