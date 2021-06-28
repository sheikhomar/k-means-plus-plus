#include <utils/random.hpp>

using namespace utils;

RandomIndexer::RandomIndexer(std::mt19937 re, size_t s) : randomEngine(re), sampler(0, s)
{
}

size_t
RandomIndexer::next()
{
    return sampler(randomEngine);
}

Random::Random(int fixedSeed)
{
    if (fixedSeed == -1)
    {
        static std::random_device randomSeed;
        this->randomEngine.seed(randomSeed());
    }
    else
    {
        this->randomEngine.seed((uint)fixedSeed);
    }
}

RandomIndexer
Random::getIndexer(size_t size)
{
    return RandomIndexer(this->randomEngine, size);
}

double
Random::getDouble()
{
    return pickRandomValue(this->randomEngine);
}

std::shared_ptr<blaze::DynamicVector<size_t>>
Random::runWeightedReservoirSampling(const size_t k, const size_t n, blaze::DynamicVector<size_t> weights)
{
    assert(weights.size() == n);

    auto indexSampler = this->getIndexer(k);
    auto data = std::make_shared<blaze::DynamicVector<size_t>>(k);

    // Algorithm by M. T. Chao
    double sum = 0;

    // Fill the reservoir array
    for (size_t i = 0; i < k; i++)
    {
        (*data)[i] = i;
        sum = sum + weights[i];
    }

    for (size_t i = k; i < n; i++)
    {
        sum = sum + weights[i];

        // Compute the probability for item i
        double p_i = k * weights[i] / sum;

        // Random value between 0 and 1
        auto q = this->getDouble();

        if (q <= p_i)
        {
            auto sampleIndex = indexSampler.next();
            (*data)[sampleIndex] = i;
        }
    }

    return data;
}
