#pragma once

#include <string>
#include <boost/array.hpp>
#include <blaze/Math.h>
#include <boost/range/algorithm_ext/erase.hpp>

namespace kmeans {
    /**
     * @brief Implementation of the k-Means clustering algorithm.
     */
    class KMeans {
        uint numOfClusters;
        bool initKMeansPlusPlus;
        uint maxIterations;
        double convergenceDiff;

    public:
        /**
         * @brief Creates a new instance of KMeans.
         * @param numOfClusters The number of clusters to generate.
         * @param initKMeansPlusPlus Initialise centroids using k-Means++.
         * @param maxIterations Maximum number of iterations.
         * @param convergenceDiff The difference in the norms of the centroids when to stop k-Means iteration.
         */
        KMeans(uint numOfClusters, bool initKMeansPlusPlus = true, uint maxIterations = 100, double convergenceDiff = 0.0001);

        /**
         * @brief Runs the algorithm.
         * @param data A NxD data matrix containing N data points of D dimensions.
         */
        blaze::DynamicVector<size_t>
        run(const blaze::DynamicMatrix<double>& data);
    };

} // namespace kmeans
