set(headers
    include/clustering/kmeans.hpp
    include/clustering/cluster_assignment_list.hpp
    include/clustering/clustering_result.hpp
    include/coresets/sensitivity_sampling.hpp
    include/utils/random.hpp
)

set(sources
    source/clustering/kmeans.cpp
    source/clustering/cluster_assignment_list.cpp
    source/clustering/clustering_result.cpp
    source/coresets/sensitivity_sampling.cpp
    source/utils/random.cpp
)

set(exe_sources
		standalone/source/main.cpp
		${sources}
)
