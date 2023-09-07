#pragma once
#define FLOAT_TYPEDEFS(T) static_assert(std::is_floating_point<T>::value, "T must be float"); typedef std::array<T,3> coord3d; typedef std::array<T,2> coord2d; typedef T real_t;
#define INT_TYPEDEFS(K) static_assert(std::is_integral<K>::value, "K must be integral type"); typedef std::array<K,3> node3; typedef std::array<K,2> node2; typedef K node_t; typedef std::array<K,6> node6;
#define TEMPLATE_TYPEDEFS(T,K) FLOAT_TYPEDEFS(T) INT_TYPEDEFS(K)
#define SEMINARIO_FORCE_CONSTANTS 0

#include "coord3d.cpp"
#include "isomer_batch.hh"
#include "sym_mat3.cpp"
#include "cubic_graph.cpp"
#include "node_neighbours.cpp"
#include "constants.cpp"
#include "matrix3.cpp"