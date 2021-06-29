#pragma once

#include <string>
#include <iostream>

#include <blaze/Math.h>
#include <boost/array.hpp>
#include <boost/range/algorithm_ext/erase.hpp>

namespace kmeans
{
    /**
     * @brief Represents a collection of points-to-cluster assignments.
     */
    class ClusterAssignmentList
    {
    public:
        /**
        * @brief Creates a new instance of ClusterAssignmentList.
        * @param numOfPoints The total number of points in the dataset.
        * @param numOfClusters The number of clusters that are generated.
        */
        ClusterAssignmentList(uint numOfPoints, uint numOfClusters);

        /**
         * @brief Assign a point to a cluster.
         * @param pointIndex The index of the point to assign the cluster to.
         * @param clusterIndex The index of the cluster to assign the point to.
         * @param distance The distance between the point and the cluster.
         */
        void
        assign(size_t pointIndex, size_t clusterIndex, double distance);

        /**
         * @brief Gets the assigned cluster index of a point.
         * @param pointIndex The index of the point for which to return the cluster index.
         */
        size_t
        getCluster(size_t pointIndex);

        /**
         * @brief Returns the total number of points in the dataset.
         */
        uint
        getNumberOfPoints();

        /**
         * @brief Returns the number of clusters that are generated.
         */
        uint
        getNumberOfClusters();

        /**
         * @brief Returns the distance of each point to its assigned cluster's centroid.
         */
        blaze::DynamicVector<double>&
        getCentroidDistances();

        /**
         * @brief Returns the number of points in a cluster.
         */
        size_t
        countPointsInCluster(size_t clusterIndex);

    private:
        /**
         * The total number of points in the dataset.
         */
        uint numOfPoints;

        /**
         * The number of clusters that are generated.
         */
        uint numOfClusters;

        /**
         * A vector of size N contain the cluster index for each point.
         */
        blaze::DynamicVector<size_t> clusters;

        /**
         * A vector of size N containing the distance between the
         * assigned cluster of each point in the dataset.
         */
        blaze::DynamicVector<double> distances;
    };

} // namespace kmeans
