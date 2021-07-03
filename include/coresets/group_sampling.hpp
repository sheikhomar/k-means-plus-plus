#pragma once

#include <algorithm>
#include <vector>
#include <iostream>
#include <cmath>

#include <clustering/kmeans.hpp>
#include <utils/random.hpp>

namespace coresets
{
    class GroupSampling
    {
    public:
        GroupSampling();
        
        void
        run(const blaze::DynamicMatrix<double> &data);
    };
}
