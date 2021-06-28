#pragma once

#include <memory>
#include <iostream>
#include <random>
#include <string>

#include <blaze/Math.h>
#include <boost/array.hpp>
#include <boost/range/algorithm_ext/erase.hpp>

#include <kmeans/cluster_assignment_list.hpp>
#include <kmeans/clustering_result.hpp>

namespace kmeans
{
    /**
     * @brief Implementation of the k-Means clustering algorithm.
     */
    class KMeans
    {
    private:
        uint numOfClusters;
        bool initKMeansPlusPlus;
        uint maxIterations;
        double convergenceDiff;
        int randomSeed;

        /**
         * @brief Creates centroids by picking points in the data matrix uniformly at random.
         * @param dataMatrix A NxD data matrix containing N data points where each point has D dimensions.
         */
        blaze::DynamicMatrix<double>
        initCentroidsNaive(const blaze::DynamicMatrix<double> &dataMatrix);

        /**
         * @brief Creates centroids according to the k-Means++ initialisation.
         * @param dataMatrix A NxD data matrix containing N data points where each point has D dimensions.
         */
        blaze::DynamicMatrix<double>
        initCentroidsKMeansPlusPlus(const blaze::DynamicMatrix<double> &dataMatrix);

        /**
         * @brief Run Lloyd's algorithm to perform the clustering of data points.
         * @param dataMatrix A NxD data matrix containing N data points where each point has D dimensions.
         * @param dataMatrix Initial k centroids where k is the number of required clusters.
         */
        std::shared_ptr<ClusteringResult>
        runLloydsAlgorithm(const blaze::DynamicMatrix<double> &dataMatrix, blaze::DynamicMatrix<double> initialCentroids);

    public:
        /**
         * @brief Creates a new instance of KMeans.
         * @param numOfClusters The number of clusters to generate.
         * @param initKMeansPlusPlus Initialise centroids using k-Means++.
         * @param maxIterations Maximum number of iterations.
         * @param convergenceDiff The difference in the norms of the centroids when to stop k-Means iteration.
         * @param randomSeed The seed for random number generators. Use a value other than -1 to fix make the algorithm deterministic.
         */
        KMeans(uint numOfClusters, bool initKMeansPlusPlus = true, uint maxIterations = 100, double convergenceDiff = 0.0001, int randomSeed = -1);

        /**
         * @brief Runs the algorithm.
         * @param data A NxD data matrix containing N data points where each point has D dimensions.
         */
        std::shared_ptr<ClusteringResult>
        run(const blaze::DynamicMatrix<double> &data);
    };

} // namespace kmeans
