#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <clustering/clustering_result.hpp>
#include <clustering/kmeans.hpp>
#include <coresets/coreset.hpp>
#include <utils/random.hpp>

namespace coresets
{
    struct InternalRing
    {
        const size_t PointIndex;
        const size_t ClusterIndex;
        const double PointCost;
        const int RangeValue;
        const double LowerBoundCost;
        const double UpperBoundCost;

        InternalRing(size_t postIndex, size_t clusterIndex, double pointCost, int rangeValue, double lowerBoundCost, double upperBoundCost) : PointIndex(postIndex), ClusterIndex(clusterIndex), PointCost(pointCost), RangeValue(rangeValue),
                                                                                                                                              LowerBoundCost(lowerBoundCost), UpperBoundCost(upperBoundCost)
        {
        }

        InternalRing &operator=(const InternalRing &) = delete; // Disallow assignment
    };

    struct ExternalRing
    {
        const size_t PointIndex;
        const size_t ClusterIndex;
        const double PointCost;
        const double CostBoundary;
        const bool IsInner;

        ExternalRing(size_t postIndex, size_t clusterIndex, double pointCost, double costBoundary, bool isInnerRing) : PointIndex(postIndex), ClusterIndex(clusterIndex), PointCost(pointCost), CostBoundary(costBoundary), IsInner(isInnerRing)
        {
        }

        ExternalRing &operator=(const ExternalRing &) = delete; // Disallow assignment
    };

    class RingSet
    {
        std::vector<std::shared_ptr<InternalRing>> internalRings;
        std::vector<std::shared_ptr<ExternalRing>> externalRings;
        const int rangeStart;
        const int rangeEnd;

    public:
        RingSet(int start, int end) : internalRings(), externalRings(), rangeStart(start), rangeEnd(end)
        {
        }

        void add(const std::shared_ptr<InternalRing> ring)
        {
            internalRings.push_back(ring);

            printf("Internal Point %3ld with cost(p, A) = %0.4f  ->  R[%2d, %ld]  [%0.4f, %0.4f) \n",
                   ring->PointIndex, ring->PointCost, ring->RangeValue, ring->ClusterIndex,
                   ring->LowerBoundCost, ring->UpperBoundCost);
        }

        void add(const std::shared_ptr<ExternalRing> ring)
        {
            externalRings.push_back(ring);

            printf("External Point %3ld with cost(p, A) = %0.4f cluster(p)=%ld -> the cost(p, A) ",
                   ring->PointIndex, ring->PointCost, ring->ClusterIndex);

            if (ring->IsInner)
            {
                printf("falls below the cost range of inner most ring (%.4f)", ring->CostBoundary);
            }
            else
            {
                printf("is above the cost range of outer most ring (%.4f)", ring->CostBoundary);
            }

            printf("\n");
        }

        size_t
        getNumberOfInnerRingPoints(size_t clusterIndex)
        {
            size_t count = 0;
            for (size_t i = 0; i < externalRings.size(); i++)
            {
                auto ring = externalRings[i];
                if (ring->ClusterIndex == clusterIndex && ring->IsInner)
                {
                    count++;
                }
            }
            return count;
        }
    };

    class GroupSampling
    {
    public:
        GroupSampling();

        void
        run(const blaze::DynamicMatrix<double> &data);

    private:
        const size_t beta;

        std::shared_ptr<RingSet>
        makeRings(const std::shared_ptr<clustering::ClusteringResult> clusters);
    };
}
