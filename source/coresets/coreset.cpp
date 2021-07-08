#include <coresets/coreset.hpp>

using namespace coresets;

Coreset::Coreset()
{
}

void Coreset::addPoint(size_t pointIndex, double weight)
{
    auto coresetPoint = std::make_shared<WeightedPoint>(pointIndex, weight, false);
    this->points.push_back(coresetPoint);
}

void Coreset::addCenter(size_t clusterIndex, double weight)
{
    auto coresetPoint = std::make_shared<WeightedPoint>(clusterIndex, weight, true);
    this->points.push_back(coresetPoint);
}

size_t
Coreset::size()
{
    return this->points.size();
}
