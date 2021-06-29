#pragma once

#include <algorithm>
#include <vector>
#include <iostream>

#include <clustering/kmeans.hpp>
#include <utils/random.hpp>

namespace coresets
{
    class CoresetResult
    {
    public:
        CoresetResult();
    };

    class SensitivySampling
    {
    public:
        SensitivySampling();
        
        std::shared_ptr<CoresetResult>
        run(const blaze::DynamicMatrix<double> &data);
    };
}
