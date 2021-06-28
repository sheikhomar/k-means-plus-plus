#pragma once

#include <kmeans/kmeans.hpp>

namespace coresets
{
    class Cora
    {
    public:
        Cora();
        void run(const blaze::DynamicMatrix<double> &data);
    };
}