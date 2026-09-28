#pragma once
#include <cplib/graph/internal/weighted_matching_engine.hpp>
namespace cplib {
// Positive edges only; parallel edges use the largest weight.
template<UnDirectedGraph G> auto maximum_weight_matching(const G& g){return detail::weighted_matching::solve<false>(g);}
}
