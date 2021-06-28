#pragma once

#include <string>
#include <iostream>
#include <random>

#include <blaze/Math.h>
#include <boost/array.hpp>
#include <boost/range/algorithm_ext/erase.hpp>

#include <kmeans/cluster_assignment_list.hpp>

namespace kmeans
{
    /**
     * @brief Represents the output of a clustering algorithm.
     */
    class ClusteringResult
    {
    public:
        /**
         * @brief Creates a new instance of ClusteringResult.
         * @param clusterAssignments Cluster assignments.
         * @param centroids The final centroids.
         */
        ClusteringResult(const kmeans::ClusterAssignmentList &clusterAssignments, blaze::DynamicMatrix<double> &centroids);

        const kmeans::ClusterAssignmentList &getClusterAssignments();

        const blaze::DynamicMatrix<double> &getCentroids();
    private:
        kmeans::ClusterAssignmentList clusterAssignments;
        blaze::DynamicMatrix<double> centroids;
    };
} // namespace kmeans
