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
        const size_t numOfDimensions;
        size_t numOfPoints;
        blaze::DynamicMatrix<double> points;
        std::vector<double> weights;

    public:
        Coreset(size_t initialSize, size_t numOfDimensions);

        void
        add(const blaze::DynamicMatrix<double> &dataPoints, size_t pointIndex, double weight);

        /**
         * @brief Returns the number of points in this coreset.
         */
        size_t
        getNumberOfPoints();
    };
}
