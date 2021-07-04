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
    public:
        std::shared_ptr<CoresetResult>
        run(const blaze::DynamicMatrix<double> &data);

        std::vector<WeightedPoint>
        calcCoresetPoints(const clustering::ClusterAssignmentList clusterAssignments, const size_t targetCoresetPoints);
    };
}
