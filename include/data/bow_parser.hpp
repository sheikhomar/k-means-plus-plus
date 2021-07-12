#pragma once

#include <algorithm>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>

#include <blaze/Math.h>

#include <boost/array.hpp>
#include <boost/range/algorithm_ext/erase.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>

namespace data
{
    /**
     * Represents a data parser.
     */
    class IDataParser
    {
    public:
        virtual ~IDataParser() {}

        /**
         * Parses the given file and converts it into a data matrix.
         */
        virtual std::shared_ptr<blaze::DynamicMatrix<double>>
        parse(const std::string &filePath) = 0; // pure virtual method
    };

    class BagOfWordsParser
    {
    public:
        std::shared_ptr<blaze::DynamicMatrix<double>>
        parse(const std::string &filePath);
    };
}
