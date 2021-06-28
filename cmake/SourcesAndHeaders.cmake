set(sources
    source/kmeans.cpp
    source/cluster_assignment_list.cpp
    source/clustering_result.cpp
)

set(exe_sources
		standalone/source/main.cpp
		${sources}
)

set(headers
    include/kmeans/kmeans.hpp
    include/kmeans/cluster_assignment_list.hpp
    include/kmeans/clustering_result.hpp
)
