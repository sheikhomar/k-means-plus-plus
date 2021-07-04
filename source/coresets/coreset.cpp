#include <coresets/coreset.hpp>

using namespace coresets;

Coreset::Coreset(size_t initialSize, size_t d) : numOfDimensions(d), numOfPoints(0), points(initialSize, d), weights()
{
}

void
Coreset::add(const blaze::DynamicMatrix<double> &dataPoints, const size_t pointIndex, const double weight)
{
    auto currentCapacity = this->points.capacity() / numOfDimensions;
    if (currentCapacity < (this->numOfPoints + 1))
    {
        auto newCapacity = ceil(static_cast<double>(currentCapacity) * 1.5);
        this->points.resize(static_cast<size_t>(newCapacity), numOfDimensions);
    }

    blaze::row(this->points, this->numOfPoints) = blaze::row(dataPoints, pointIndex);
    weights.push_back(weight);
    numOfPoints++;
}

size_t
Coreset::getNumberOfPoints()
{
    return this->numOfPoints;
}
