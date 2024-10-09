#pragma once
#include <fullerenes/graph.hh>
#include <fullerenes/polyhedron.hh>
#include <fullerenes/sycl-headers/sycl-status-enum.hh>
#include <fullerenes/sycl-headers/sycl-fullerene-misc-tuple-fun.hh>
#include <fullerenes/sycl-headers/reference-wrapper.hh>


template <typename T = float, typename K = uint16_t>
struct Fullerene
{
    static_assert(std::is_floating_point<T>::value, "T must be a floating point type");
    static_assert(std::is_integral<K>::value, "K must be an integral type");

    Fullerene(const FullereneDataMembers<Span, T, K>& data, const FullereneMetaMembers<ReferenceWrapper, K>& meta, size_t N, size_t Nf) 
        : d_(data), m_(meta), N_(N), Nf_(Nf) {}

    explicit operator Graph() const {
        ConditionFunctor cond(StatusFlag::FULLERENEGRAPH_PREPARED);
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
        ConditionFunctor cubic_and_3d(StatusFlag::FULLERENEGRAPH_PREPARED | StatusFlag::CONVERGED_3D);
        bool is_polyhedron = (int)(m_.flags_.get()) & (int)(StatusFlag::FULLERENEGRAPH_PREPARED);
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

    Fullerene<T, K> &operator=(const Graph &G) {
        if (this->N_ == 0 || this->Nf_ == 0) {throw std::invalid_argument("Fullerenes are non-owning, cannot assign to an uninitialized fullerene");}
        if (G.neighbours.size() != N_ && G.neighbours.size() != Nf_) {throw std::invalid_argument("Graph has incompatible number of vertices: " + std::to_string(G.neighbours.size()) + " vs N: " + std::to_string(N_) + " or Nf: " + std::to_string(Nf_));}
        bool is_cubic = G.neighbours.size() == N_ && G.neighbours[0].size() == 3;
        auto A = is_cubic ? d_.A_cubic_.template as_span<K>() : d_.A_dual_.template as_span<K>();
        auto count = is_cubic ? N_ : Nf_;
        for (size_t i = 0; i < count; i++) {
            auto degree = is_cubic ? 3 : G.neighbours[i].size();
            if(!is_cubic) {d_.deg_[i] = degree;}
            for (size_t j = 0; j < degree; j++) {
                A[i * (is_cubic ? 3 : 6) + j] = G.neighbours[i][j];
            }
        }
        m_.flags_.get() |= (is_cubic ? StatusFlag::FULLERENEGRAPH_PREPARED : StatusFlag::DUAL_INITIALIZED);
        return *this;
    }

    Fullerene<T, K>& operator=(const std::tuple<std::reference_wrapper<const Graph>, size_t>& Graph_and_ID){
        auto& [G, ID] = Graph_and_ID;
        *this = G.get();
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