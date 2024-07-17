#pragma once
#include <string>
#include <array>
#include <span>
#include <memory>
#include "graph.hh"
#include "polyhedron.hh"


//Use this definition in class definitions to add the following data members

enum class Policy
{
    SYNC,
    ASYNC
};

enum class DeviceType
{
    CPU,
    GPU,
    ACCELERATOR,
    HOST,
    NUM_DEV_TYPES
};

enum class DeviceProperty
{
    MAX_WORK_GROUP_SIZE,
    MAX_CLOCK_FREQUENCY,
    MAX_COMPUTE_UNITS,
    MAX_MEM_ALLOC_SIZE,
    GLOBAL_MEM_SIZE, 
    LOCAL_MEM_SIZE,
    MAX_CONSTANT_ARGS,
    NUMBER_OF_PROPERTIES
};
typedef size_t device_id_t;
typedef std::pair<device_id_t, DeviceType> Device;

namespace fullerene_detail {
    template<typename... Args1, typename... Args2, std::size_t ... Is>
    auto merge_tuples(std::tuple<Args1&...> t1, std::tuple<Args2&...> t2,
                      std::index_sequence<Is...>) {
        return std::forward_as_tuple(std::tie(std::get<Is>(t1), std::get<Is>(t2))...);
    }

    template<typename Tuple, std::size_t... Is>
    auto rm_const_helper(Tuple&& t, std::index_sequence<Is...>) {
        return std::forward_as_tuple(const_cast<std::remove_const_t<std::tuple_element_t<Is, std::remove_reference_t<Tuple>>>&>(std::get<Is>(std::forward<Tuple>(t)))...);
    }

    template<typename... Args>
    auto rm_const(std::tuple<Args...>&& t) {
        return rm_const_helper(std::move(t), std::index_sequence_for<Args...>{});
    }

    template <std::size_t ...I, typename F> auto with_sequence_impl(F &&func, std::index_sequence<I...>)
    {
        return func(std::integral_constant<std::size_t, I>{}...);
    }

    template <std::size_t N, typename F> auto with_sequence(F &&func)
    {
        return with_sequence_impl(std::forward<F>(func), std::make_index_sequence<N>{});
    }
}

template<typename... Args1, typename... Args2>
auto forward_merge_tuples(const std::tuple<Args1&...>& t1, const std::tuple<Args2&...>& t2) {
    static_assert(sizeof...(Args1) == sizeof...(Args2));
    return fullerene_detail::merge_tuples(t1, t2, std::make_index_sequence<sizeof...(Args1)>());
}

// SyclContext Must be created after any fork() calls
struct SyclContext
{
    SyclContext();
    ~SyclContext();

    size_t device_get_property(Device dev, DeviceProperty property) const;
    size_t device_get_count(DeviceType device_type) const;
    size_t device_get_count() const;
    std::string device_get_name(Device dev) const;

    bool queue_is_in_order(size_t queue_id) const;
    void queue_wait(size_t queue_id) const;

    void queue_push_back(Device dev, bool in_order = true);
    size_t queue_get_count() const;

private:
    const std::vector<std::byte> cpus_;
    const std::vector<std::byte> gpus_;
    const std::vector<std::byte> accelerators_;
    std::vector<std::byte> queues_;
};

struct SyclQueueImpl;

struct SyclQueue{

    /*  Must be declared here and defined in the implementation file otherwise the 
        compiler will attempt to generate a default constructor, 
        which involves default constructing the unique_ptr, 
        which in turn requires the definition of SyclQueueImpl*/
    SyclQueue(); 
    SyclQueue(Device device, bool in_order = true);
    SyclQueue(std::unique_ptr<SyclQueueImpl>&& impl);
    ~SyclQueue();

    
    
    void wait() const;
    void wait_and_throw() const;
    
    std::unique_ptr<SyclQueueImpl> impl_;
    
    const Device device_;
    const bool in_order_;
};

enum class StatusFlag
{
    EMPTY = 1 << 0,            // 0000 0001
    CONVERGED_2D = 1 << 1,     // 0000 0010
    CONVERGED_3D = 1 << 2,     // 0000 0100
    PLZ_CHECK = 1 << 3,        // 0000 1000
    FAILED_2D = 1 << 4,        // 0001 0000
    FAILED_3D = 1 << 5,        // 0010 0000
    DUAL_INITIALIZED = 1 << 6, // 0100 0000
    CUBIC_INITIALIZED = 1 << 7 // 1000 0000
};

template <class T>
inline T operator~(T a) { return (T) ~(int)a; }
template <class T, class K>
inline T operator|(T a, K b) { return (T)((int)a | (int)b); }
template <class T, class K>
inline T operator&(T a, K b) { return (T)((int)a & (int)b); }
template <class T, class K>
inline T operator^(T a, K b) { return (T)((int)a ^ (int)b); }
template <class T, class K>
inline T &operator|=(T &a, K b) { return (T &)((int &)a |= (int)b); }
template <class T, class K>
inline T &operator&=(T &a, K b) { return (T &)((int &)a &= (int)b); }
template <class T, class K>
inline T &operator^=(T &a, K b) { return (T &)((int &)a ^= (int)b); }

template <class T, class K>
bool not_set(const T flag, const K condition) { return int(flag & condition) == 0; }
template <class T, class K>
bool all_set(const T flag, const K condition) { return int(flag | ~condition) == ~0; }
template <class T, class K>
bool any_set(const T flag, const K condition) { return int(flag & condition) != 0; }

struct ConditionFunctor
{
    ConditionFunctor() : not_conditions(0), and_conditions(0), or_conditions(~0) {}                                                                                                                                                                    // Default, no conditions
    ConditionFunctor(int and_conditions, int not_conditions = 0, int or_conditions = ~0) : not_conditions(not_conditions), and_conditions(and_conditions), or_conditions(or_conditions) {}                                                             // Custom conditions;
    ConditionFunctor(StatusFlag and_conditions, StatusFlag not_conditions = (StatusFlag)0, StatusFlag or_conditions = (StatusFlag)~0) : not_conditions((int)not_conditions), and_conditions((int)and_conditions), or_conditions((int)or_conditions) {} // Custom conditions;
    
    ConditionFunctor(const ConditionFunctor &other) = default;
    ConditionFunctor(ConditionFunctor &&other) = default;
    ConditionFunctor &operator=(const ConditionFunctor &other) = default;
    ConditionFunctor &operator=(ConditionFunctor &&other) = default;

    const int not_conditions;                                                                                                                                                                                                                          // Bitwise condtions that must not be met
    const int and_conditions;                                                                                                                                                                                                                          // Bitwise condtions that must be met
    const int or_conditions;                                                                                                                                                                                                                           // Bitwise condtions of which at least one must be met

    // First checks if the flags that must be met are met, then checks if the flags that must not be met are not met
    inline constexpr bool operator()(const StatusFlag flag) const
    {
        int iflag = (int)flag;
        return not_set(iflag, not_conditions) && all_set(iflag, and_conditions) && any_set(iflag, or_conditions);
    }
    inline constexpr bool operator()(const int flag) const { return not_set(flag, not_conditions) && all_set(flag, and_conditions) && any_set(flag, or_conditions); }
};
template <typename T>
struct SyclVectorIterator
{
    explicit SyclVectorIterator(T *ptr) : data_(ptr) {}
    T &operator*() const { return *data_; }
    T *operator->() const { return data_; }
    SyclVectorIterator &operator++()
    {
        ++data_;
        return *this;
    }
    SyclVectorIterator operator++(int)
    {
        SyclVectorIterator tmp = *this;
        ++data_;
        return tmp;
    }
    SyclVectorIterator &operator--()
    {
        --data_;
        return *this;
    }
    SyclVectorIterator operator--(int)
    {
        SyclVectorIterator tmp = *this;
        --data_;
        return tmp;
    }

    friend bool operator==(const SyclVectorIterator &a, const SyclVectorIterator &b) { return a.data_ == b.data_; };
    friend bool operator!=(const SyclVectorIterator &a, const SyclVectorIterator &b) { return a.data_ != b.data_; };
    friend size_t operator-(const SyclVectorIterator &a, const SyclVectorIterator &b) { return a.data_ - b.data_; };
    friend bool operator<(const SyclVectorIterator &a, const SyclVectorIterator &b) { return a.data_ < b.data_; };
    friend bool operator>(const SyclVectorIterator &a, const SyclVectorIterator &b) { return a.data_ > b.data_; };

private:
    T *data_;
};

template <typename T>
struct Span
{
    inline constexpr Span() : data_(nullptr), size_(0) {}
    inline constexpr Span(T *data, size_t size) : data_(data), size_(size) {}
    inline constexpr Span(T *begin, T *end) : data_(begin), size_(std::distance(begin, end)) {}
    inline constexpr Span(const Span<T> &other) = default;
    inline constexpr Span(Span<T> &&other) = default;
    inline constexpr Span<T> subspan(size_t offset, size_t count) const { return Span<T>(data_ + offset, count); }
    inline constexpr Span<T>& operator= (const Span<T> &other) { data_ = other.data_; size_ = other.size_; return *this; }
    inline constexpr Span<T>& operator= (Span<T> &&other) { return *this = other; }
    inline bool operator==(const Span<T> &other) const;
    inline constexpr T &operator[](size_t index) const { return data_[index]; }
    inline constexpr T &at(size_t index) const{assert(index < size_); return data_[index];}
    inline constexpr T *data() const { return data_; }
    inline constexpr size_t size() const { return size_; }
    inline constexpr bool empty() const { return size_ == 0; }
    inline constexpr size_t size_bytes() const { return size_ * sizeof(T); }
    inline constexpr T *begin() const { return data_; }
    inline constexpr T *end() const { return data_ + size_; }
    inline constexpr T &front() const { return data_[0]; }
    inline constexpr T &back() const { return data_[size_ - 1]; }
    
private:
    T *data_;
    size_t size_;
};

namespace fullerene_detail{
    template<typename... Args1, typename... Args2, std::size_t ... Is, size_t N>
    auto construct_spans_impl(size_t idx, std::array<int, N>&& size_factors, std::tuple<Args1&...> dst, std::tuple<Args2&...> src,
                        std::index_sequence<Is...>) {
        auto assign_span = [](auto& lhs, auto rhs_ptr, size_t size) {
            lhs = Span(rhs_ptr, size);
        };

        (..., assign_span(std::get<Is>(dst),std::get<Is>(src).data() + size_factors[Is] * ((int)idx), size_factors[Is]));
    }

    template<size_t N, typename... Args1, typename... Args2>
    void construct_spans(size_t idx, std::array<int, N>&& size_factors, std::tuple<Args1&...> dst, std::tuple<Args2&...> src) {
        construct_spans_impl(idx, std::move(size_factors), dst, src, std::make_index_sequence<sizeof...(Args1)>{});
    }
}



template <typename T>
struct SyclVector
{
    SyclVector();
    SyclVector(size_t size);
    SyclVector(size_t size, T value);
    SyclVector(const SyclVector<T> &other);
    SyclVector(SyclVector<T> &&other);
    SyclVector<T> &operator=(const SyclVector<T> &other);
    SyclVector<T> &operator=(SyclVector<T> &&other);

    ~SyclVector();

    operator Span<T>() const;
    Span<T> operator()(size_t offset, size_t count) const;

    void fill(T data);
    T *data() const;
    size_t size() const;
    size_t capacity() const;
    void resize(size_t new_size);
    void reserve(size_t new_capacity);
    void clear();

    T &operator[](size_t index);
    const T &operator[](size_t index) const;

    T &at(size_t index);
    const T &at(size_t index) const;

    bool operator==(const SyclVector<T> &other) const;

    void push_back(const T &value);
    void pop_back();

    template <typename U>
    friend std::ostream &operator<<(std::ostream &os, const SyclVector<U> &vec);

    T *begin() const;
    T *end() const;

    inline constexpr T &back() const { return data_[size_ - 1]; }
    inline constexpr T &front() const { return data_[0]; }

private:
    size_t size_;
    size_t capacity_;
    T *data_;
};


template <template<typename> class Container, typename T, typename K>
struct FullereneDataMembers{
    mutable Container<T> X_cubic_; // 3D Embedding of the Cubic Graph {N * 3}
    mutable Container<T> X_dual_;  // 3D Embedding of the Dual Graph {Nf * 6}
    mutable Container<K> A_cubic_; // Adjacency Matrix (Cubic) {N * 3}
    mutable Container<K> A_dual_;  // Adjacency Matrix (Dual) {Nf * 6}
    mutable Container<K> faces_;   // Atom indices of the faces {Nf * 6}
    mutable Container<K> deg_;     // Vertex degrees in the dual graph, face degrees in the cubic graph, face degrees in the dual graph is always 3 (it is a triangulation)
    // Container<K> quad_edge_;   //Quad edge representation of the mapping between the cubic and dual graphs {Nf * 6}

    FullereneDataMembers() = default;
    ~FullereneDataMembers() = default;
    FullereneDataMembers<Container, T, K>(const FullereneDataMembers<Container, T, K> &other) = default;
    FullereneDataMembers<Container, T, K>(FullereneDataMembers<Container, T, K> &&other) = default;
    FullereneDataMembers<Container, T, K> &operator=(const FullereneDataMembers<Container, T, K> &other) = default;
    FullereneDataMembers<Container, T, K> &operator=(FullereneDataMembers<Container, T, K> &&other) = default;

    inline constexpr auto to_tuple() const { return std::forward_as_tuple(X_cubic_, X_dual_, A_cubic_, A_dual_, faces_, deg_); }
    
    static inline constexpr auto get_size_factors(int N, int capacity) { 
            int Nf = N/2 + 2;
            return std::array{(int)N*3*capacity, (int)Nf*3*capacity, (int)N*3*capacity, (int)Nf*6*capacity, (int)Nf*6*capacity, (int)Nf*capacity};
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

    
    inline constexpr auto to_tuple() const { return std::forward_as_tuple(ID_, iterations_, flags_, valid_indices_); }
    static inline constexpr auto get_size_factors(int capacity) { return std::array{(int)capacity, (int)capacity, (int)capacity, (int)capacity}; }
};

template <typename T = float, typename K = uint16_t>
struct Fullerene
{
    static_assert(std::is_floating_point<T>::value, "T must be a floating point type");
    static_assert(std::is_integral<K>::value, "K must be an integral type");

    Fullerene(const FullereneDataMembers<Span, T, K>& data, size_t N, size_t Nf);

    Fullerene(const Graph &G, bool is_cubic = false);
    Fullerene(const neighbours_t &neighbours, bool is_cubic = false);
    Fullerene(const Polyhedron &P);

    explicit operator Graph() const;
    explicit operator Polyhedron() const;

    bool operator==(const Fullerene &other) const;
    ~Fullerene() = default;

    FullereneDataMembers<Span, T, K> d_;

    const size_t N_;                      // Number of vertices in the cubic graph {1}
    const size_t Nf_;                     // Number of faces in the dual graph {1}
};

template <typename T = float, typename K = uint16_t>
struct FullereneBatch
{
    static_assert(std::is_floating_point<T>::value, "T must be a floating point type");
    static_assert(std::is_integral<K>::value, "K must be an integral type");

    FullereneBatch();
    ~FullereneBatch() = default;

    FullereneBatch(size_t N, int capacity);

    void push_back(const Graph &G, const int ID = -1);
    void push_back(const neighbours_t &neighbours, bool is_cubic = false, const int ID = -1);
    void push_back(const Polyhedron &P, const int ID = -1);
    void push_back(const Fullerene<T, K> &fullerene, const int ID = -1);
    void resize(size_t new_capacity);

    FullereneBatch(const FullereneBatch<T, K> &other) = default;
    FullereneBatch(FullereneBatch<T, K> &&other) = default;
    FullereneBatch<T, K> &operator=(const FullereneBatch<T, K> &other) = default;

    Fullerene<T, K> operator[](size_t index) const;

    // Getters for private members that are not supposed to be modified outside the class
    constexpr inline int size()     const {return size_;}
    constexpr inline int capacity() const {return capacity_;}
    constexpr inline int front()    const {return front_;}
    constexpr inline int back()     const {return back_;}

    FullereneDataMembers<SyclVector, T, K> d_;
    FullereneMetaMembers<SyclVector, K> m_;

    const size_t N_;                      // Number of vertices in the cubic graph {1}
    const size_t Nf_;                     // Number of faces in the dual graph {1}

protected:
    int capacity_; // Maximum number of isomers in the batch {1}
    int size_;     // Number of isomers in the batch {1}
    int front_;    // Index of the first isomer in the batch {1}
    int back_;     // Index of the last isomer in the batch {1}

    void prepare_for_push_back(const neighbours_t &neighbours, bool is_cubic);
};

// FullereneQueue and FullereneBatch are functionally the we use the type system to express intent and to prevent misuse,
// FullereneQueues have the additional property that their non-StatusFlag::EMPTY isomers are contiguously stored in memory
template <typename T = float, typename K = uint16_t>
struct FullereneQueue : public FullereneBatch<T, K>
{
    FullereneQueue() = default;
    ~FullereneQueue() = default;

    FullereneQueue(size_t N, int capacity) : FullereneBatch<T, K>(N, capacity) {}
    FullereneQueue(const FullereneQueue<T, K> &other) = default;
    FullereneQueue(FullereneQueue<T, K> &&other) = default;
    FullereneQueue<T, K>& operator=(const FullereneQueue<T, K> &other) = default;
    FullereneQueue<T, K>& operator=(FullereneQueue<T, K> &&other) = default;
    
    inline static ConditionFunctor empty_cond = ConditionFunctor(StatusFlag::EMPTY, ~StatusFlag::EMPTY); //Checks that the isomer is EMPTY and that it is not any other status
    inline static ConditionFunctor initialized_cond = ConditionFunctor(0,0,(int) (StatusFlag::DUAL_INITIALIZED | StatusFlag::CUBIC_INITIALIZED)); //Checks that the isomer is either DUAL_INITIALIZED or CUBIC_INITIALIZED
    inline static ConditionFunctor converged_cond = ConditionFunctor(StatusFlag::CONVERGED_2D | StatusFlag::CONVERGED_3D); //Checks that the isomer is both CONVERGED_2D and CONVERGED_3D

    friend void push(SyclQueue& Q, FullereneQueue<T, K> &dst_queue, FullereneBatch <T, K>&src_batch, ConditionFunctor transfer_cond);
    friend void push(SyclQueue& Q, FullereneBatch<T, K> &dst_batch, FullereneQueue<T, K> &src_queue, ConditionFunctor transfer_cond);

    private: 
        //This will be used internally to implement the push functions
        friend void push_impl(SyclQueue& Q, FullereneBatch<T, K> &dst_batch, FullereneBatch<T, K> &src_batch, ConditionFunctor transfer_cond);
};

template <typename T = float, typename K = uint16_t>
struct FullereneBatchView
{
    static_assert(std::is_floating_point<T>::value, "T must be a floating point type");
    static_assert(std::is_integral<K>::value, "K must be an integral type");

    FullereneBatchView() = default;
    ~FullereneBatchView() = default;

    FullereneBatchView(const FullereneBatchView<T, K> &other) = default;
    FullereneBatchView(FullereneBatchView<T, K> &&other) = default;

    FullereneBatchView(const FullereneBatch<T, K> &batch, size_t offset = 0, int count = -1) : 
        N_(batch.N_), Nf_(batch.Nf_),
        size_(count < 0 ? batch.size() - offset : count){
        if (offset < 0 || (offset + count) >= batch.size()) {throw std::out_of_range("Offset out of range");}
        fullerene_detail::construct_spans(offset, d_.get_size_factors(batch.N_, size_), d_.to_tuple(), batch.d_.to_tuple());
        fullerene_detail::construct_spans(offset, m_.get_size_factors(size_), m_.to_tuple(), batch.m_.to_tuple());
    } 

    Fullerene<T, K> operator[](size_t index) const {
        if (index >= size_) {assert(!"Index out of range");}
        return Fullerene<T, K>(d_, N_, Nf_);
    }
                                
    // Getters for private members that are not supposed to be modified outside the class
    constexpr inline int size() const {return size_;}                            // Returns the size of the view
    constexpr inline Fullerene<T, K> front() const { return *this[0]; }         // Returns the first fullerene in the view
    constexpr inline Fullerene<T, K> back() const { return *this[size() - 1]; } // Returns the last fullerene in the view
    
    FullereneDataMembers<Span, T, K> d_;
    FullereneMetaMembers<Span, K> m_;

    const size_t N_;                      // Number of vertices in the cubic graph {1}
    const size_t Nf_;                     // Number of faces in the dual graph {1}
    const int size_; // Number of isomers in the batch {1}
};
struct Identity{
    template <typename T>
    constexpr T&& operator()(T&& x) const noexcept { return std::forward<T>(x);}
};

struct Square{
    template <typename T>
    constexpr T operator()(T x) const noexcept { return x*x;}
};

struct Plus{
    template <typename T>
    constexpr T operator()(T x, T y) const noexcept { return x + y;}
};

struct Minus{
    template <typename T>
    constexpr T operator()(T x, T y) const noexcept { return x - y;}
};

#define COMMA ,

#define DECLARE_FUNCS(templates, returntype, funcname, sig_args, call_args) \
    templates returntype funcname(SyclQueue& Q, T* begin, T* end, T* store, sig_args);\
    templates returntype funcname(SyclQueue& Q, const SyclVector<T> &vec, SyclVector<T> &result, sig_args);\
    templates returntype funcname(SyclQueue& Q, const Span<T> vec, const Span<T> result, sig_args);\

#define DECLARE_OR_DEFINE_ALL_FUNCS(definemode)\
        definemode(template <typename T COMMA typename BinaryOp>, void, exclusive_scan, T init COMMA BinaryOp op, init COMMA op)\
        definemode(template <typename T COMMA typename BinaryOp>, void, inclusive_scan, T init COMMA BinaryOp op, init COMMA op)\
        definemode(template <typename T COMMA typename BinaryOp COMMA typename UnaryOp>, void, transform_exclusive_scan, T init COMMA BinaryOp op COMMA UnaryOp f, init COMMA op COMMA f)\
        definemode(template <typename T COMMA typename BinaryOp COMMA typename UnaryOp>, void, transform_inclusive_scan, T init COMMA BinaryOp op COMMA UnaryOp f, init COMMA op COMMA f)\
        definemode(template <typename T COMMA typename UnaryOp>, void, transform, UnaryOp f, f)\
        definemode(template <typename T COMMA typename BinaryOp>, T, reduce, T init COMMA BinaryOp op, init COMMA op)\
        definemode(template <typename T COMMA typename BinaryOp COMMA typename UnaryOp>, T, transform_reduce, T init COMMA BinaryOp op COMMA UnaryOp f, init COMMA op COMMA f)\
        //definemode(template <typename T>, void, fill, Policy policy = Policy::SYNC)\
        //definemode(template <typename T>, void, unique, Policy policy = Policy::SYNC)\
        //definemode(template <typename T>, void, sort, Policy policy = Policy::SYNC)\
        //definemode(template <typename T>, void, reverse, Policy policy = Policy::SYNC)
namespace primitives{
    DECLARE_OR_DEFINE_ALL_FUNCS(DECLARE_FUNCS)
}


