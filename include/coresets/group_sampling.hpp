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

        ClusteredPoint &operator=(const ClusteredPoint &) = delete; // Disallow assignment
    };

    /**
     * Represents a group which is uniquely identified by its range value (j) and its ring range value (l) i.e., G_{j,l}
     */
    class Group
    {

    public:
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

        Group(size_t rangeValue, int ringRangeValue, double lowerBoundCost, double upperBoundCost) : RangeValue(rangeValue), RingRangeValue(ringRangeValue), LowerBoundCost(lowerBoundCost), UpperBoundCost(upperBoundCost)
        {
        }

        Group &operator=(const Group &) = delete; // Disallow assignment

        void addPoint(size_t point, size_t cluster, double cost)
        {
            points.push_back(std::make_shared<ClusteredPoint>(point, cluster, cost));
        }

        const std::vector<std::shared_ptr<ClusteredPoint>> &
        getPoints() const
        {
            return points;
        }

        double
        calcTotalCost() const
        {
            double sum = 0;
            for (size_t i = 0; i < points.size(); i++)
            {
                sum += points[i]->Cost;
            }

            return sum;
        }

    private:
        /**
         * The points assigned to this group.
         */
        std::vector<std::shared_ptr<ClusteredPoint>> points;
    };

    class GroupSet
    {
        std::vector<std::shared_ptr<Group>> groups;

    public:
        GroupSet &operator=(const GroupSet &) = delete; // Disallow assignment

        std::shared_ptr<Group> create(size_t rangeValue, int ringRangeValue, double lowerBoundCost, double upperBoundCost)
        {
            auto group = std::make_shared<Group>(rangeValue, ringRangeValue, lowerBoundCost, upperBoundCost);
            groups.push_back(group);
            return group;
        }

        size_t size() const
        {
            return this->groups.size();
        }

        std::shared_ptr<Group>
        operator[](size_t index) const
        {
            return this->groups[index];
        }

        std::shared_ptr<Group>
        get(size_t index) const
        {
            return this->groups[index];
        }

        blaze::DynamicVector<double>
        calcNormalizedCosts() const
        {
            blaze::DynamicVector<double> costs(this->groups.size());
            double sumOfGroupCosts = 0.0;
            for (size_t i = 0; i < this->groups.size(); i++)
            {
                auto groupCost = this->groups[i]->calcTotalCost();
                costs[i] = groupCost;
                sumOfGroupCosts += groupCost;
            }

            return costs / sumOfGroupCosts;
        }
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

    class Ring
    {
    public:
        /**
         * The cluster for which this ring belongs to.
         */
        const size_t ClusterIndex;

        /**
         * The range value of this ring i.e., `l`.
         */
        const int RangeValue;

        /**
         * The average cost for the cluster associated with this ring.
         */
        const double AverageClusterCost;

        Ring(size_t clusterIndex, int rangeValue, double averageClusterCost) : ClusterIndex(clusterIndex), RangeValue(rangeValue), AverageClusterCost(averageClusterCost)
        {
            // Ring upper bound cost := Δ_c * 2^l
            LowerBoundCost = averageClusterCost * std::pow(2, rangeValue);

            // Ring upper bound cost := Δ_c * 2^(l+1)
            UpperBoundCost = averageClusterCost * std::pow(2, rangeValue + 1);

            TotalCost = 0.0;
        }

        Ring &operator=(const Ring &) = delete; // Disallow assignment

        /**
         * @brief Adds a point to the ring if its costs is within the bounds of this ring.
         * @param pointIndex The index of the point to add to this ring.
         * @param cost The cost of the point i.e., cost(p, A)
         * @return `true` if points is added to the ring, otherwise `false`
         */
        bool
        tryAddPoint(size_t pointIndex, double cost)
        {
            if (isCostWithinBounds(cost))
            {
                printf("Internal Point %3ld with cost(p, A) = %0.4f  ->  R[%2d, %ld]  [%0.4f, %0.4f) \n",
                       pointIndex, cost, RangeValue, ClusterIndex, LowerBoundCost, UpperBoundCost);

                points.push_back(std::make_shared<ClusteredPoint>(pointIndex, ClusterIndex, cost));
                TotalCost += cost;
                return true;
            }

            return false;
        }

        bool
        isCostWithinBounds(double cost)
        {
            // If cost(p, A) is between Δ_c*2^l and Δ_c*2^(l+1) ...
            return cost >= LowerBoundCost && cost < UpperBoundCost;
        }

        double getLowerBoundCost() { return LowerBoundCost; }
        double getUpperBoundCost() { return UpperBoundCost; }

        /**
         * @brief Sums the costs of all points in captured by this ring.
         */
        double getTotalCost() { return TotalCost; }

        const std::vector<std::shared_ptr<ClusteredPoint>> &
        getPoints() const
        {
            return this->points;
        }

    private:
        /**
         * The points assigned to this ring.
         */
        std::vector<std::shared_ptr<ClusteredPoint>> points;

        double LowerBoundCost;

        double UpperBoundCost;

        /**
         * The sum of the point cost in this ring.
         */
        double TotalCost;
    };

    class RingSet
    {
        std::vector<std::shared_ptr<ExternalRing>> externalRings;
        std::vector<std::shared_ptr<Ring>> rings;

    public:
        const int RangeStart;
        const int RangeEnd;

        RingSet(int start, int end) : externalRings(), RangeStart(start), RangeEnd(end)
        {
        }

        RingSet &operator=(const RingSet &) = delete; // Disallow assignment

        std::shared_ptr<Ring> find(size_t clusterIndex, int rangeValue) const
        {
            for (size_t i = 0; i < rings.size(); i++)
            {
                auto ring = rings[i];
                if (ring->ClusterIndex == clusterIndex && ring->RangeValue == rangeValue)
                {
                    return ring;
                }
            }
            return nullptr;
        }

        std::shared_ptr<Ring> findOrCreate(size_t clusterIndex, int rangeValue, double averageClusterCost)
        {
            auto ring = find(clusterIndex, rangeValue);
            if (ring == nullptr)
            {
                // printf("Ring for cluster=%ld and l=%2d not found. Creating...\n", clusterIndex, rangeValue);
                ring = std::make_shared<Ring>(clusterIndex, rangeValue, averageClusterCost);
                rings.push_back(ring);
            }
            return ring;
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

        /**
         * @brief Sums the costs of all points in captured by all rings for a given range i.e., cost(R_l) = sum_{p in R_l} cost(p, A)
         * @param ringRangeValue The ring range value i.e. l
         */
        double
        calcRingCost(int ringRangeValue) const
        {
            double sum = 0.0F;
            for (size_t i = 0; i < rings.size(); i++)
            {
                auto ring = rings[i];
                if (ring->RangeValue == ringRangeValue)
                {
                    sum += ring->getTotalCost();
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
        makeRings(const clustering::ClusterAssignmentList &clusters);

        /**
         * @brief Add points inside doughnut holes i.e., points that are closest to cluster centers but are not captured by any rings.
         */
        void addShortfallPointsToCoreset(const clustering::ClusterAssignmentList &clusters, const std::shared_ptr<RingSet> rings, std::vector<WeightedPoint> &coresetPoints);

        /**
         * @brief Add points outside doughnuts i.e., points that are far from cluster centers and are not captured by any rings.
         */
        void addOvershootPointsToCoreset(const blaze::DynamicMatrix<double> &data, std::shared_ptr<RingSet> rings, std::vector<WeightedPoint> &coresetPoints);

        std::shared_ptr<GroupSet>
        makeGroups(const clustering::ClusterAssignmentList &clusters, const std::shared_ptr<RingSet> rings, const size_t numberOfGroups);
    };
}
