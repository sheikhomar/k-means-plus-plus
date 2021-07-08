#pragma once

#include <algorithm>
#include <vector>
#include <iostream>

#include <clustering/kmeans.hpp>
#include <utils/random.hpp>

namespace coresets
{
    struct WeightedPoint
    {
        const size_t Index;
        const double Weight;
        const bool IsCenter;

        WeightedPoint(size_t index, double weight, bool isCenter) : Index(index), Weight(weight), IsCenter(isCenter)
        {
        }
    };

    class Coreset
    {
        std::vector<std::shared_ptr<WeightedPoint>> points;

    public:
        Coreset();

        /**
         * @brief Adds a point to the coreset.
         * @param pointIndex The index of the point to add to the coreset.
         * @param weight The weight of the point.
         */
        void addPoint(size_t pointIndex, double weight);

        /**
         * @brief Adds a center to the coreset.
         * @param clusterIndex The index of the center to add to the coreset.
         * @param weight The weight of the center.
         */
        void addCenter(size_t clusterIndex, double weight);

        /**
         * @brief Returns the number of points in this coreset.
         */
        size_t
        size();
    };
}
