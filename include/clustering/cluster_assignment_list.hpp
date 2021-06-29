#pragma once

#include <string>
#include <iostream>

#include <blaze/Math.h>
#include <boost/array.hpp>
#include <boost/range/algorithm_ext/erase.hpp>

namespace clustering
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
        * @brief Copy constructor.
        * @param other The other cluster assignments to copy from.
        */
        ClusterAssignmentList(const ClusterAssignmentList &other);

        /**
         * @brief Copies cluster assignments from another object.
         * @param other The other cluster assignments to copy from.
         */
        ClusterAssignmentList&
        operator=(const ClusterAssignmentList &other);

        /**
         * @brief Assign all data points to their closest centers.
         */
        void
        assignAll(const blaze::DynamicMatrix<double> &dataPoints, const blaze::DynamicMatrix<double> &centers);

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
        blaze::DynamicVector<double> &
        getCentroidDistances();

        /**
         * @brief Returns the number of points in a cluster.
         */
        size_t
        countPointsInCluster(size_t clusterIndex);

        /**
         * @brief Returns the cost of the cluster assignments. 
         * 
         * The cost is the sum of pairwise distances between each point and its closest center.
         */
        double
        calcCost();

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

}
