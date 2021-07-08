#include <coresets/coreset.hpp>

using namespace coresets;

Coreset::Coreset(size_t targetSize) : TargetSize(targetSize)
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

std::shared_ptr<WeightedPoint>
Coreset::at(size_t index) const
{
    return this->points[index];
}

size_t
Coreset::size() const
{
    return this->points.size();
}
