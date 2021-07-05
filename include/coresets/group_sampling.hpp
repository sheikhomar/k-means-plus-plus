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
    /**
     * Represents a reference to a point which has been assigned a cluster.
     */
    struct ClusteredPoint
    {
        /**
         * The index of the point in the dataset.
         */
        const size_t PointIndex;

        /**
         * The index of the cluster for which this point is assigned.
         */
        const size_t ClusterIndex;

        /**
         * The cost of this point in its assigned cluster. 
         * 
         * The cost is the distance to the assigned cluster's center.
         */
        const double Cost;

        ClusteredPoint(size_t pointIndex, size_t clusterIndex, double cost) : PointIndex(pointIndex), ClusterIndex(clusterIndex), Cost(cost)
        {
        }
    };

    /**
     * Represents a group which is uniquely identified by its range value (j) and its ring range value (l) i.e., G_{j,l}
     */
    class Group
    {

    public:
        Group(size_t rangeValue, int ringRangeValue, double lowerBoundCost, double upperBoundCost) : RangeValue(rangeValue), RingRangeValue(ringRangeValue), LowerBoundCost(lowerBoundCost), UpperBoundCost(upperBoundCost)
        {
        }

        void addPoint(size_t point, size_t cluster, double cost)
        {
            points.push_back(std::make_shared<ClusteredPoint>(point, cluster, cost));
        }

    private:
        /**
         * The group's range value `j`, it is a non-negative value.
         */
        const size_t RangeValue;

        /**
         * The group's ring range value i.e., `l`.
         */
        const int RingRangeValue;

        /**
         * The lower bound cost of the group.
         */
        const double LowerBoundCost;

        /**
         * The upper bound cost of the group.
         */
        const double UpperBoundCost;

        /**
         * The points assigned to this group.
         */
        std::vector<std::shared_ptr<ClusteredPoint>> points;
    };

    class GroupSet
    {
        std::vector<std::shared_ptr<Group>> groups;

    public:
        std::shared_ptr<Group> create(size_t rangeValue, int ringRangeValue, double lowerBoundCost, double upperBoundCost)
        {
            auto group = std::make_shared<Group>(rangeValue, ringRangeValue, lowerBoundCost, upperBoundCost);
            groups.push_back(group);
            return group;
        }
    };

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

    public:
        const int RangeStart;
        const int RangeEnd;

        RingSet(int start, int end) : internalRings(), externalRings(), RangeStart(start), RangeEnd(end)
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
        countInternalRings() const
        {
            return this->internalRings.size();
        }

        /**
         * @brief Sums the costs of all points in captured by ring for a given range i.e., cost(R_l) = sum_{p in R_l} cost(p, A)
         * @param ringRangeValue The ring range value i.e. l
         */
        double
        calcRingCost(int ringRangeValue) const
        {
            double sum = 0.0F;
            for (size_t i = 0; i < internalRings.size(); i++)
            {
                auto ring = internalRings[i];
                if (ring->RangeValue == ringRangeValue)
                {
                    sum += ring->PointCost;
                }
            }
            return sum;
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

        std::vector<size_t>
        getPointsOutsideAllRings() const
        {
            std::vector<size_t> pointsOutsideAllRings;
            for (size_t i = 0; i < externalRings.size(); i++)
            {
                auto ring = externalRings[i];
                if (ring->IsInner == false)
                {
                    pointsOutsideAllRings.push_back(ring->PointIndex);
                }
            }

            return pointsOutsideAllRings;
        }

        std::vector<std::shared_ptr<InternalRing>>
        getInternalRingPointsInCluster(size_t clusterIndex, int ringRangeValue)
        {
            std::vector<std::shared_ptr<InternalRing>> rings;

            for (size_t i = 0; i < this->internalRings.size(); i++)
            {
                auto ring = internalRings[i];
                if (ring->ClusterIndex == clusterIndex && ring->RangeValue == ringRangeValue)
                {
                    rings.push_back(ring);
                }
            }

            return rings;
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

        void addInnerMostRingPoints(const clustering::ClusterAssignmentList &clusters, const std::shared_ptr<RingSet> rings, std::vector<WeightedPoint> &coresetPoints);

        void addPointsOutsideAllRings(const blaze::DynamicMatrix<double> &data, std::shared_ptr<RingSet> rings, std::vector<WeightedPoint> &coresetPoints);

        std::shared_ptr<GroupSet>
        makeGroups(const clustering::ClusterAssignmentList &clusters, const std::shared_ptr<RingSet> rings, const size_t numberOfGroups);
    };
}
