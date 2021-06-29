set(headers
    include/clustering/cluster_assignment_list.hpp
    include/clustering/clustering_result.hpp
    include/clustering/local_search.hpp
    include/clustering/kmeans.hpp
    include/coresets/sensitivity_sampling.hpp
    include/utils/random.hpp
)

set(sources
    source/clustering/cluster_assignment_list.cpp
    source/clustering/clustering_result.cpp
    source/clustering/local_search.cpp
    source/clustering/kmeans.cpp
    source/coresets/sensitivity_sampling.cpp
    source/utils/random.cpp
)

set(exe_sources
		standalone/source/main.cpp
		${sources}
)
