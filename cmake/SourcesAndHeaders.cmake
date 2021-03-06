set(headers
    include/clustering/cluster_assignment_list.hpp
    include/clustering/clustering_result.hpp
    include/clustering/local_search.hpp
    include/clustering/kmeans.hpp
    include/coresets/coreset.hpp
    include/coresets/group_sampling.hpp
    include/coresets/sensitivity_sampling.hpp
    include/coresets/stream_km.hpp
    include/data/bow_parser.hpp
    include/data/census_parser.hpp
    include/data/covertype_parser.hpp
    include/data/data_parser.hpp
    include/data/tower_parser.hpp
    include/utils/random.hpp
)

set(sources
    source/clustering/cluster_assignment_list.cpp
    source/clustering/clustering_result.cpp
    source/clustering/local_search.cpp
    source/clustering/kmeans.cpp
    source/coresets/coreset.cpp
    source/coresets/group_sampling.cpp
    source/coresets/sensitivity_sampling.cpp
    source/coresets/stream_km.cpp
    source/data/bow_parser.cpp
    source/data/census_parser.cpp
    source/data/covertype_parser.cpp
    source/data/tower_parser.cpp
    source/utils/random.cpp
)

set(exe_sources
		standalone/source/main.cpp
		${sources}
)
