#pragma once

#include <algorithm>
#include <iostream>
#include <memory>
#include <random>
#include <vector>

#include <blaze/Math.h>
#include <boost/array.hpp>

namespace utils
{
    class RandomIndexer
    {
    public:
        RandomIndexer(std::mt19937 randomEngine, size_t size);
        size_t next();

    private:
        std::mt19937 randomEngine;
        std::uniform_int_distribution<size_t> sampler;
    };

    class Random
    {
    public:
        Random(int fixedSeed);

        /**
         * Returns a random real number in the interval [0.0, 1.0).
         */
        double
        getDouble();

        RandomIndexer
        getIndexer(size_t size);

        std::shared_ptr<blaze::DynamicVector<size_t>>
        runWeightedReservoirSampling(const size_t k, const size_t n, blaze::DynamicVector<size_t> weights);

        /**
         * @brief Randomly select `k` indices from an array of size `n` with replacement.
         * @param k The number of indices to pick from.
         * @param weight A collection of weights associated with each entry.
         */
        std::shared_ptr<blaze::DynamicVector<size_t>>
        choice(const size_t k, const size_t n, blaze::DynamicVector<size_t> weights);

        /**
         * @brief Randomly select an index using the given weights.
         */
        size_t
        choice(blaze::DynamicVector<size_t> weights);

    private:
        std::mt19937 randomEngine;
        std::uniform_real_distribution<> pickRandomValue;
    };
}
