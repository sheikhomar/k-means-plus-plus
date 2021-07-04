#pragma once

#include <algorithm>
#include <vector>
#include <iostream>

#include <clustering/kmeans.hpp>
#include <coresets/coreset.hpp>
#include <utils/random.hpp>

namespace coresets
{
    class CoresetResult
    {
    public:
        CoresetResult();
    };

    class SensitivitySampling
    {
        utils::Random random;

    public:
        SensitivitySampling(const int randomSeed);

        std::shared_ptr<CoresetResult>
        run(const blaze::DynamicMatrix<double> &data);

        std::vector<WeightedPoint>
        calcCoresetPoints(const std::shared_ptr<clustering::ClusteringResult> result, const size_t targetCoresetPoints);
    };
}
