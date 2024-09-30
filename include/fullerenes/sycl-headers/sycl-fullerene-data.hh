#pragma once
#include <tuple>
#include <functional>
#include <array>

template <template<typename> class Container, typename T, typename K>
struct FullereneDataMembers{
    mutable Container<T> X_cubic_; // 3D Embedding of the Cubic Graph {N * 3}
    mutable Container<T> X_dual_;  // 3D Embedding of the Dual Graph {Nf * 6}
    mutable Container<K> A_cubic_; // Adjacency Matrix (Cubic) {N * 3}
    mutable Container<K> A_dual_;  // Adjacency Matrix (Dual) {Nf * 6}
    mutable Container<K> faces_cubic_;   // Atom indices of the hexagons/pentagons {Nf * 6}
    mutable Container<K> faces_dual_;    // "Face-indices" of the triangles {N * 3}
    mutable Container<K> deg_;     // Vertex degrees in the dual graph, face degrees in the cubic graph, face degrees in the dual graph is always 3 (it is a triangulation)

    FullereneDataMembers() = default;
    ~FullereneDataMembers() = default;
    FullereneDataMembers<Container, T, K>(const FullereneDataMembers<Container, T, K> &other) = default;
    FullereneDataMembers<Container, T, K>(FullereneDataMembers<Container, T, K> &&other) = default;
    FullereneDataMembers<Container, T, K> &operator=(const FullereneDataMembers<Container, T, K> &other) = default;
    FullereneDataMembers<Container, T, K> &operator=(FullereneDataMembers<Container, T, K> &&other) = default;
    bool operator==(const FullereneDataMembers<Container, T, K> &other) const {
        auto compute_equality = [](auto... args) { return ((args.first == args.second) && ...); };
        return compute_equality(std::make_pair(to_tuple(), other.to_tuple()));
    }

    inline constexpr auto to_tuple() const 
        { return std::make_tuple(std::ref(X_cubic_), 
            std::ref(X_dual_), 
            std::ref(A_cubic_), 
            std::ref(A_dual_), 
            std::ref(faces_cubic_), 
            std::ref(faces_dual_),
            std::ref(deg_)); }
    
    static inline constexpr auto get_size_factors(int N, int capacity) { 
            int Nf = N/2 + 2;
            return std::array{(int)N*3*capacity,    //X_cubic_ 
                            (int)Nf*3*capacity,     //X_dual_
                            (int)N*3*capacity,      //A_cubic_
                            (int)Nf*6*capacity,     //A_dual_
                            (int)Nf*6*capacity,     //faces_cubic_
                            (int)N*3*capacity,      //faces_dual_
                            (int)Nf*capacity};      //deg_
    }
};

template <template<typename> class Container, typename K>
struct FullereneMetaMembers{
    mutable Container<size_t> ID_;               // Buckygen ID of the isomer {1}
    mutable Container<K> iterations_;       // Number of forcefield CG iterations performed so far {1}
    mutable Container<StatusFlag> flags_;   // Status flags of the isomers {1}
    mutable Container<K> valid_indices_;    // Indices of the valid isomers {1}

    FullereneMetaMembers() = default;
    ~FullereneMetaMembers() = default;
    FullereneMetaMembers<Container, K>(const FullereneMetaMembers<Container, K> &other) = default;
    FullereneMetaMembers<Container, K>(FullereneMetaMembers<Container, K> &&other) = default;
    FullereneMetaMembers<Container, K> &operator=(const FullereneMetaMembers<Container, K> &other) = default;
    FullereneMetaMembers<Container, K> &operator=(FullereneMetaMembers<Container, K> &&other) = default;

    bool operator==(const FullereneMetaMembers<Container, K> &other) const {
        auto compute_equality = [](auto... args) { return ((args.first == args.second) && ...); };
        return compute_equality(std::make_pair(to_tuple(), other.to_tuple()));
    }
    
    inline constexpr auto to_tuple() const { return std::make_tuple(std::ref(ID_), std::ref(iterations_), std::ref(flags_), std::ref(valid_indices_)); }
    static inline constexpr auto get_size_factors(int capacity) { return std::array{(int)capacity, (int)capacity, (int)capacity, (int)capacity}; }
};