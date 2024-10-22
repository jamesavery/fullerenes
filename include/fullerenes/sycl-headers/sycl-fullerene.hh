#pragma once
#include <fullerenes/graph.hh>
#include <fullerenes/polyhedron.hh>
#include <fullerenes/sycl-headers/sycl-status-enum.hh>
#include <fullerenes/sycl-headers/sycl-fullerene-misc-tuple-fun.hh>
#include <fullerenes/sycl-headers/sycl-fullerene-data.hh>
#include <fullerenes/sycl-headers/reference-wrapper.hh>

namespace extra_type_traits {
    template <typename T>
    struct is_array_t : std::false_type {};

    template <typename T>
    struct is_array_t<std::vector<T>> : std::true_type {};

    template <typename T>
    struct is_array_t<Span<T>> : std::true_type {};

    template <typename T, size_t N>
    struct is_array_t<std::array<T, N>> : std::true_type {};

    template <typename T>
    struct is_tuple_t : std::false_type {};

    template <typename... Ts>
    struct is_tuple_t<std::tuple<Ts...>> : std::true_type {};

    template <typename T>
    struct is_pair_t : std::false_type {};

    template <typename T, typename U>
    struct is_pair_t<std::pair<T, U>> : std::true_type {};

    // Generalization to handle references and const qualifiers
    template <typename T>
    struct is_array_t<T&> : is_array_t<T> {};

    template <typename T>
    struct is_array_t<T&&> : is_array_t<T> {};

    template <typename T>
    struct is_array_t<const T> : is_array_t<T> {};
}

using namespace condition_detail;
template <typename T = float, typename K = uint16_t>
struct Fullerene
{   
    static_assert(std::is_floating_point<T>::value, "T must be a floating point type");
    static_assert(std::is_integral<K>::value, "K must be an integral type");

    Fullerene(const FullereneDataMembers<Span, T, K>& data, const FullereneMetaMembers<ReferenceWrapper, K>& meta, size_t N, size_t Nf) 
        : d_(data), m_(meta), N_(N), Nf_(Nf) {}

    explicit operator Graph() const {
        ConditionFunctor cond(StatusEnum::FULLERENEGRAPH_PREPARED);
        bool is_cubic = cond(m_.flags_.get());
        Graph G(neighbours_t(is_cubic ? N_ : Nf_));
        auto A = is_cubic ? d_.A_cubic_.template as_span<K>() : d_.A_dual_.template as_span<K>();
        auto count = is_cubic ? N_ : Nf_;
        for (size_t i = 0; i < count; i++) {
            auto degree = is_cubic ? 3 : d_.deg_[i];
            for (size_t j = 0; j < degree; j++) {
                G.neighbours[i].push_back(A[i * (is_cubic ? 3 : 6) + j]);
            }
        }
        return G;
    }
    
    explicit operator Polyhedron() const {
        using namespace condition_detail;
        if (N_ == 0 || Nf_ == 0) {throw std::invalid_argument("Fullerenes are non-owning, cannot convert to Polyhedron without initializing the fullerene");}
        const ConditionFunctor graph_condition(0, 0, StatusEnum::FULLERENEGRAPH_PREPARED | StatusEnum::DUAL_INITIALIZED | StatusEnum::CUBIC_INITIALIZED);
        const ConditionFunctor geometry_condition(0, 0, StatusEnum::NOT_CONVERGED | StatusEnum::CONVERGED_3D | StatusEnum::FAILED_3D);
        bool is_polyhedron = graph_condition(m_.flags_.get()) && geometry_condition(m_.flags_.get());
        if (!is_polyhedron) {throw std::invalid_argument("Fullerene is not a valid polyhedron, Flag: " + std::to_string(m_.flags_.get()));}
        Polyhedron P;
        using points_t = decltype(P.points);
        points_t points(N_);
        neighbours_t neighbours(N_);
        auto& A = d_.A_cubic_;
        for (size_t i = 0; i < N_; i++) {
            for (size_t j = 0; j < 3; j++) {
                neighbours[i].push_back(A[i][j]);
                points[i][j] = d_.X_cubic_[i][j];
            }
        }
        return Polyhedron(PlanarGraph(Graph(neighbours,true)), points);
    }




    bool operator==(const Fullerene other) const;

    Fullerene(const Fullerene<T, K> &other) = default;
    Fullerene<T, K> &operator=(const Fullerene<T, K> &other) = default;

    Fullerene<T, K> &operator=(const neighbours_t &neighbours) {
        if (this->N_ == 0 || this->Nf_ == 0) {throw std::invalid_argument("Fullerenes are non-owning, cannot assign to an uninitialized fullerene");}
        if (neighbours.size() != N_ && neighbours.size() != Nf_) {throw std::invalid_argument("Graph has incompatible number of vertices: " + std::to_string(neighbours.size()) + " vs N: " + std::to_string(N_) + " or Nf: " + std::to_string(Nf_));}
        bool is_cubic = neighbours.size() == N_ && neighbours[0].size() == 3;
        auto A = is_cubic ? d_.A_cubic_.template as_span<K>() : d_.A_dual_.template as_span<K>();
        auto count = is_cubic ? N_ : Nf_;
        for (size_t i = 0; i < count; i++) {
            auto degree = is_cubic ? 3 : neighbours[i].size();
            if(!is_cubic) {d_.deg_[i] = degree;}
            for (size_t j = 0; j < degree; j++) {
                A[i * (is_cubic ? 3 : 6) + j] = neighbours[i][j];
            }
        }
        m_.flags_.get() = (is_cubic ? StatusEnum::CUBIC_INITIALIZED : StatusEnum::DUAL_INITIALIZED);
        m_.iterations_.get() = 0;
        return *this;
    }

    Fullerene<T, K> &operator=(std::tuple<std::reference_wrapper<const neighbours_t>, std::reference_wrapper<const vector<coord2d>>> neighbours_and_layout) {
        auto &[neighbours, layout] = neighbours_and_layout;
        *this = neighbours.get();
        auto is_cubic = N_ == neighbours.get().size();
        auto dst_ptr = is_cubic ? d_.X_cubic_.template as_span<std::array<T, 2>>() : d_.X_dual_.template as_span<std::array<T, 2>>();
        //using LayoutType = std::decay_t<decltype(layout.get()[0])>;
        
        auto copy_layout = [](auto& dst, auto& src) {
            static_assert(extra_type_traits::is_array_t<decltype(src)>::value, "Source layout: expected an array type");
            static_assert(extra_type_traits::is_array_t<decltype(dst)>::value, "Destination layout: expected an array type");
            std::size(src) == std::size(dst) ? std::size(src) : throw std::invalid_argument("Layouts have different sizes");

            if constexpr (extra_type_traits::is_array_t<decltype(dst[0])>::value && extra_type_traits::is_array_t<decltype(src[0])>::value) {
                if constexpr (std::is_same_v<decltype(dst[0][0]), decltype(src[0][0])>) {
                    memcpy(dst.data(), src.data(), src.size() * sizeof(src[0])); // Fast path for same types.
                } else {
                    for (size_t i = 0; i < src.size(); i++) {
                        for (size_t j = 0; j < 2; j++) {
                            dst[i][j] = src[i][j];
                        }
                    }
                }
            } else if constexpr (extra_type_traits::is_pair_t<decltype(src[0])>::value && extra_type_traits::is_array_t<decltype(dst[0])>::value) {
                for (size_t i = 0; i < src.size(); i++) {
                    dst[i][0] = static_cast<T>(src[i].first);
                    dst[i][1] = static_cast<T>(src[i].second);
                }
            } else if constexpr (extra_type_traits::is_pair_t<decltype(src[0])>::value && extra_type_traits::is_pair_t<decltype(dst[0])>::value) {
                for (size_t i = 0; i < src.size(); i++) {
                    dst[i].first = static_cast<T>(src[i].first);
                    dst[i].second = static_cast<T>(src[i].second);
                }
            } else {
                throw std::invalid_argument("Unsupported layout types");
                
            }
        };
        copy_layout(dst_ptr, layout.get());
        m_.flags_.get() |= StatusEnum::CONVERGED_2D;
        return *this;
    }

    Fullerene<T, K> &operator=(const Graph &G) {
        *this = G.neighbours;
        return *this;
    }

    Fullerene<T, K> &operator=(const PlanarGraph &PG) {
        *this = std::make_tuple(std::cref(PG.neighbours), std::cref(PG.layout2d));
        return *this;
    }

    Fullerene<T, K> &operator=(const FullereneGraph &FG) {
        *this = std::make_tuple(std::cref(FG.neighbours), std::cref(FG.layout2d));
        return *this;
    }

    Fullerene<T, K> &operator=(const FullereneDual &FD) {
        *this = std::make_tuple(std::cref(FD.neighbours), std::cref(FD.layout2d));
        return *this;
    }

    Fullerene<T, K> &operator=(const Polyhedron &P) {
        *this = P.neighbours;
        auto is_cubic = P.points.size() == N_;
        auto dst_ptr = is_cubic ? d_.X_cubic_.data() : d_.X_dual_.data();
        if constexpr (std::is_same_v<decltype(d_.X_cubic_[0]), decltype(P.points[0])>) {
            memcpy(dst_ptr, P.points.data(), P.points.size() * sizeof(P.points[0])); //Fast path for same types.
        } else {
            for (size_t i = 0; i < P.points.size(); i++) {
                for (size_t j = 0; j < 3; j++) {
                    dst_ptr[i][j] = P.points[i][j];
                }
            }
        }
        m_.flags_.get() |= StatusEnum::CONVERGED_3D; //TODO: Ought we to check this first?
        return *this;
    }

    Fullerene<T, K> &operator=(const std::tuple<std::reference_wrapper<const neighbours_t>, size_t> &neighbours_and_ID) {
        auto &[neighbours, ID] = neighbours_and_ID;
        *this = neighbours.get();
        m_.ID_.get() = ID;
        return *this;
    }

    Fullerene<T, K>& operator=(const std::tuple<std::reference_wrapper<const Graph>, size_t>& Graph_and_ID){
        auto& [G, ID] = Graph_and_ID;
        *this = G.get().neighbours;
        m_.ID_.get() = ID;
        return *this;
    }

    Fullerene<T, K>& operator=(const std::tuple<std::reference_wrapper<const PlanarGraph>, size_t>& PlanarGraph_and_ID){
        auto& [PG, ID] = PlanarGraph_and_ID;
        *this = PG.get().neighbours;
        m_.ID_.get() = ID;
        return *this;
    }

    Fullerene<T, K>& operator=(const std::tuple<std::reference_wrapper<const FullereneGraph>, size_t>& FullereneGraph_and_ID){
        auto& [FG, ID] = FullereneGraph_and_ID;
        *this = FG.get().neighbours;
        m_.ID_.get() = ID;
        return *this;
    }

    Fullerene<T, K>& operator=(const std::tuple<std::reference_wrapper<const Polyhedron>, size_t>& Polyhedron_and_ID){
        auto& [P, ID] = Polyhedron_and_ID;
        *this = P.get();
        m_.ID_.get() = ID;
        return *this;
    }

    ~Fullerene() = default;

    template <typename U, typename V>
    static inline void copy(Fullerene<U, V> dst, const Fullerene<T, K> src){
        if (dst.N_ != src.N_ || dst.Nf_ != src.Nf_) {assert(!"Fullerenes have different sizes");}

        auto tuple_pair = forward_merge_tuples(dst.d_.to_tuple(), src.d_.to_tuple());
        
        auto copy_pair = [](auto& dst_range, auto& src_range) {
            std::copy(std::begin(src_range), std::end(src_range), std::begin(dst_range));
        };

        std::apply([&](auto&... args) { 
            (copy_pair(std::get<0>(args), std::get<1>(args)), ...); 
        }, tuple_pair);

        auto tuple_pair_meta = forward_merge_tuples(dst.m_.to_tuple(), src.m_.to_tuple());
        auto copy_meta_pair = [](auto& pair) {
            std::get<0>(pair).get() = std::get<1>(pair).get();
        };
        std::apply([&](auto&... args) { 
            (copy_meta_pair(args), ...); 
        }, tuple_pair_meta);
    }

    const FullereneDataMembers<Span, T, K> d_;
    const FullereneMetaMembers<ReferenceWrapper, K> m_;

    const size_t N_ = 0;                      // Number of vertices in the cubic graph {1}
    const size_t Nf_ = 0;                     // Number of faces in the dual graph {1}
    //const ReferenceWrapper<StatusFlag> flag_; //Reference Wrapper to the status flag of the isomer {1}
                                               //Safer than using a pointer, and is trivially copyable as opposed to a reference
                                               //Trivial copyability is required in SYCL kernel contexts

    //Output stream operator for Fullerene
    template <typename U, typename V>
    friend std::ostream& operator<<(std::ostream& os, const Fullerene<U, V>& fullerene);
};