#include <fullerenes/sycl-headers/sycl-parallel-primitives.hh>
#include <fullerenes/sycl-headers/reference-wrapper.hh>
#include <fullerenes/sycl-headers/sycl-fullerene.hh>
#include <fullerenes/sycl-headers/sycl-fullerene-structs.hh>
#include "primitives.cc"
#include <CL/sycl.hpp>
#include <unordered_map>
#include <exception>
#include <limits>
#include <cstdint>
#include "coord3d.cc"

#ifndef DEVICE_CAST
    #define DEVICE_CAST(x,ix) (reinterpret_cast<const sycl::device*>(x)[ix])
#endif
#include "queue-impl.cc"
//#define HOST_ALLOCATOR_QUEUE reinterpret_cast<sycl::queue*>(QueueWrapper::allocator_queue_.data())[0]

using namespace sycl;

template <typename U>
constexpr std::ostream& operator<<(std::ostream& os, const ReferenceWrapper<U>& ref) {
    os << ref.get();
    return os;
}

template <typename U, size_t N>
constexpr std::ostream& operator<<(std::ostream& os, const std::array<U,N>& arr) {
    os << "[";
    for (size_t i = 0; i < N; i++) {
        os << arr[i];
        if(i < N - 1) os << ", ";
    }
    os << "]";
    return os;
}

template <typename T>
SyclVector<T>::SyclVector(size_t size) : size_(size), capacity_(size) {
    data_ = sycl::malloc_shared<T>(size, sycl::device(default_selector_v), sycl::context(device(default_selector_v)));
}

template <typename T>
SyclVector<T>::SyclVector(size_t size, T value) : size_(size), capacity_(size) {
    data_ = sycl::malloc_shared<T>(size, sycl::device(default_selector_v), sycl::context(device(default_selector_v)));
    for(size_t i = 0; i < size; i++) data_[i] = value;
}

template <typename T>
SyclVector<T>::SyclVector(const SyclVector<T>& other) : size_(other.size_), capacity_(other.capacity_) {
    data_ = sycl::malloc_shared<T>(capacity_, sycl::device(default_selector_v), sycl::context(device(default_selector_v))); 
    for(size_t i = 0; i < size_; i++) data_[i] = other.data_[i];
}

template <typename T>
SyclVector<T>::~SyclVector() {if(data_) sycl::free(data_,  sycl::context(device(default_selector_v)));}

template <typename T>
void SyclVector<T>::resize(size_t new_size) {
    if(new_size > capacity_){
        T* new_data = sycl::malloc_shared<T>(new_size, sycl::device(default_selector_v), sycl::context(device(default_selector_v)));
        memcpy(new_data, data_, size_*sizeof(T));
        if(data_) sycl::free(data_,  sycl::context(device(default_selector_v)));
        data_ = new_data;
        capacity_ = new_size;
    }
    size_ = new_size;
}

template <typename T>
void SyclVector<T>::resize(size_t new_size, T val){
    if(new_size > capacity_){
        T* new_data = sycl::malloc_shared<T>(new_size, sycl::device(default_selector_v), sycl::context(device(default_selector_v)));
        memcpy(new_data, data_, size_*sizeof(T));
        if(data_) sycl::free(data_,  sycl::context(device(default_selector_v)));
        data_ = new_data;
        capacity_ = new_size;
        std::generate(data_ + size_ , data_ + new_size, [val](){return val;});
        size_ = new_size;
    }
}

template <typename T>
void SyclVector<T>::resize(size_t new_size, size_t front, size_t back, size_t seg_size){
    if(new_size > capacity_){
        T* new_data = sycl::malloc_shared<T>(new_size, sycl::device(default_selector_v), sycl::context(device(default_selector_v)));
        memset(new_data, 0, new_size*sizeof(T));
        if (capacity_ > 0){
            auto n_first_segment = back < front ? capacity_ - front : (back - front + seg_size);
            assert(n_first_segment <= new_size);
            memcpy(new_data, data_ + front, n_first_segment*sizeof(T));
            assert((n_first_segment + (back + 1)) <= new_size);
            if( back < front) memcpy(new_data + n_first_segment, data_, (back + seg_size)*sizeof(T));
        }
        if(data_) sycl::free(data_,  sycl::context(device(default_selector_v)));
        data_ = new_data;
        size_ = new_size;
        capacity_ = new_size;
    }
}

template <typename T>
void SyclVector<T>::reserve(size_t new_capacity) {
    if(new_capacity > capacity_){
        T* new_data = sycl::malloc_shared<T>(new_capacity, sycl::device(default_selector_v), sycl::context(device(default_selector_v)));
        memcpy(new_data, data_, size_*sizeof(T));
        if(data_) sycl::free(data_,  sycl::context(device(default_selector_v)));
        data_ = new_data;
        capacity_ = new_capacity;
    }
}

template <typename T>
SyclVector<T>& SyclVector<T>::operator=(const SyclVector<T>& other) {
    if(data_) sycl::free(data_,  sycl::context(device(default_selector_v)));
    size_ = other.size_;
    capacity_ = other.capacity_;
    data_ = sycl::malloc_shared<T>(capacity_, sycl::device(default_selector_v), sycl::context(device(default_selector_v)));
    for(size_t i = 0; i < size_; i++) data_[i] = other.data_[i];
    return *this;
}

template <typename U>
std::ostream& operator<<(std::ostream& os, const SyclVector<U>& vec) {
    os << (Span<U>)vec;
    return os;
}

template <typename U>
std::ostream& operator<<(std::ostream& os, const Span<U>& vec) {
    os << "[" ;
     for (size_t i = 0; i < vec.size(); i++) {
        os << vec[i];
        if(i < vec.size() - 1) os << ", ";
    }
    os << "]";
    return os;
}

template <typename T>
bool Span<T>::operator==(const Span<T> other) const {
    if(size_ != other.size_) return false;    
    if(data() == other.data()) return true;
    if constexpr (is_floating_point_v<T>){
        return std::equal(begin(), end(), other.begin(), [](auto& a,auto& b){
            T eps = std::numeric_limits<T>::epsilon() * 10;
            T max_v = std::max<T>(std::abs<T>(a), std::abs<T>(b));
            return std::abs<T>(a - b) / (max_v > eps ? max_v : 1) < eps;});
    } else{
        return std::equal(begin(), end(), other.begin());
    }
    return true;
}



template <typename T, typename K>
bool Fullerene<T,K>::operator==(const Fullerene<T,K> other) const {
    return (d_ == other.d_ && m_ == other.m_&& N_ == other.N_ && Nf_ == other.Nf_);
}

template <typename U, typename V>
std::ostream& operator<<(std::ostream& os, const Fullerene<U, V>& fullerene) {
    os << "Fullerene with " << fullerene.N_ << " vertices and " << fullerene.Nf_ << " faces\n";
    os << "Data Members:\n";
    auto print_data = [&os](auto... args) { (..., (os << args << "\n")); };
    std::apply(print_data, fullerene.d_.to_tuple());
    os  << "ID: " << fullerene.m_.ID_ << "\tFlag: " << fullerene.m_.flags_ << "\tIterations: " << fullerene.m_.iterations_ << "\n";

    return os;
}

template <typename U, typename V>
std::ostream& operator<<(std::ostream& os, const FullereneBatchView<U, V>& batch) {
    os << "FullereneBatchView with " << batch.size_ << " isomers\n";
    os << "N: " << batch.N_ << " Nf: " << batch.Nf_ << "\n";
    os << "Meta Data:\n";
    os << "Data Members:\n";
    std::for_each(batch.begin(), batch.end(), [&os, j = 0](auto fullerene) mutable { os << "Batch[ " << j++ << " ]\n" << fullerene; });
    return os;
}

template <typename U, typename V>
std::ostream& operator<<(std::ostream& os, const FullereneBatch<U, V>& batch) { os << FullereneBatchView<U, V>(batch); return os; }

template <typename T, typename K>
FullereneBatch<T,K>::FullereneBatch() : N_(0), Nf_(0), capacity_(0), size_(0) {}


template <size_t N, typename... Args, std::size_t... I>
void resize_all(std::array<int,N>&& sizes, std::tuple<Args...>&& args, std::index_sequence<I...>){
    (std::get<I>(args).resize(sizes[I], typename std::decay_t<decltype(std::get<I>(args))>::value_type{}), ...);
}

template <typename T, typename K>
FullereneBatch<T,K>::FullereneBatch(size_t N, int capacity) :
    N_(N),
    Nf_(N/2 + 2),
    capacity_(capacity),
    size_(0)
{
    resize_all(d_.get_size_factors(N_, capacity_), d_.to_tuple(), std::make_index_sequence<std::tuple_size_v<decltype(d_.to_tuple())>>{});
    resize_all(m_.get_size_factors(capacity_), m_.to_tuple(), std::make_index_sequence<std::tuple_size_v<decltype(m_.to_tuple())>>{});
}

template <typename T, typename K>
void FullereneBatch<T,K>::resize(size_t new_capacity) {
    if(new_capacity <= capacity_) return;
    capacity_ = new_capacity;
    resize_all(d_.get_size_factors(N_, capacity_), d_.to_tuple(), std::make_index_sequence<std::tuple_size_v<decltype(d_.to_tuple())>>{});
    resize_all(m_.get_size_factors(capacity_), m_.to_tuple(), std::make_index_sequence<std::tuple_size_v<decltype(m_.to_tuple())>>{});
}

template <typename T, typename K>
void FullereneQueue<T,K>::resize(size_t new_capacity) {
    if(new_capacity <= capacity_) return;
    capacity_ = new_capacity;
    auto resize_all_circularly = [&]<std::size_t... I>(std::index_sequence<I...>, auto&& tup, auto&& size_factors){
        (std::get<I>(tup).resize(size_factors[I]*new_capacity, size_factors[I]*front_, size_factors[I]*back_, size_factors[I]), ...);
    };
    resize_all_circularly(std::make_index_sequence<std::tuple_size_v<decltype(d_.to_tuple())>>{}, d_.to_tuple(), d_.get_size_factors(N_, 1));
    resize_all_circularly(std::make_index_sequence<std::tuple_size_v<decltype(m_.to_tuple())>>{}, m_.to_tuple(), m_.get_size_factors(1));       
    front_ = 0;
    back_ = size_ - 1;
}
template <typename T, typename K>
void FullereneBatch<T,K>::prepare_for_push_back(const neighbours_t& neighbours, bool is_cubic){
    if(N_ == 0 || Nf_ == 0) {
        if(is_cubic) {N_ = neighbours.size(); Nf_ = N_/2 + 2;}
        else{ Nf_ = neighbours.size(); N_ = 2*Nf_ - 4;}
    } else if((is_cubic && N_ != neighbours.size()) || (!is_cubic && Nf_ != neighbours.size())) throw std::runtime_error("Graph size is incompatible with the batch");
    if(capacity_ == size_) resize(capacity_ > 0 ? 2*capacity_ : 1);
}

template <typename T, typename K>
void FullereneBatch<T,K>::push_back(const neighbours_t& neighbours, bool is_cubic, const int ID) {
    prepare_for_push_back(neighbours, is_cubic);
    size_t N = is_cubic ? N_ : Nf_;
    for(size_t i = 0; i < N; i++){
        if(!is_cubic) d_.deg_[i + size_*Nf_] = neighbours[i].size();
        for(size_t j = 0; j < neighbours[i].size(); j++){
            if(is_cubic){d_.A_cubic_[i + size_ * N][j] = neighbours[i][j];}
            else{d_.A_dual_[i + size_*Nf_][j] = neighbours[i][j];}
        }
    }
    m_.ID_[size_] = ID;
    m_.flags_[size_] = is_cubic ? StatusEnum::FULLERENEGRAPH_PREPARED : StatusEnum::DUAL_INITIALIZED;
    m_.iterations_[size_] = 0;
    size_++;
}

template <typename T, typename K>
void FullereneBatch<T,K>::push_back(const Graph& G, const int ID) {push_back(G.neighbours, false, ID);}

template <typename T, typename K>
void FullereneBatch<T,K>::push_back(const PlanarGraph& G, const int ID) {push_back(G.neighbours, true, ID);}

template <typename T, typename K>
void FullereneBatch<T,K>::push_back(const Polyhedron& P, const int ID) {
    push_back(P.neighbours, true, ID);
    //If floating point types are not the same, we need to convert
    if(!std::is_same<T, decltype(P.points[0][0])>::value){
        
    for(size_t i = 0; i < N_; i++){
        for(size_t j = 0; j < 3; j++){
            d_.X_cubic_[i + (size_-1)*N_][j] = P.points[i][j];
        }
    } //Otherwise a contiguous memory copy will do
    } else{ memcpy(d_.X_cubic_.data() + (size_-1)*N_, P.points.data(), N_*sizeof(std::array<T,3>));}
    m_.flags_[(size_-1)] |= StatusEnum::CONVERGED_3D;
    std::cout << "Pushed Polyhedron with ID: " << ID << std::endl;
    std::cout << "Pushed Polyhedron with flags: " << m_.flags_[(size_-1)] << std::endl;
}

//Parameter pack expanded copy function
template <typename... Args>
void copy_spans(sycl::group<1>& grp, Args&&... args) {
    auto bdim = grp.get_local_linear_range();
    auto tix = grp.get_local_linear_id();
    ([&] {
        for (int i = tix; i < std::size(std::get<0>(args)); i += bdim) std::get<1>(args).data()[i] = std::get<0>(args).data()[i];
    }(), ...);
}

template <typename T, typename K> struct FullereneBatchFill{};
//Use enum class to express source and destination conditions like:
//fill(... , !GeometryStatus::EMPTY)

template <typename T> struct is_fullerene_queue : std::false_type{};
template <typename T, typename K> struct is_fullerene_queue<FullereneQueue<T,K>> : std::true_type{};

template <typename T, typename K>
void QueueUtil::push(SyclQueue& Q, FullereneQueue<T, K> &dst_queue, FullereneBatch <T, K>&src_batch, ConditionFunctor transfer_cond){
    push_impl(Q, dst_queue, src_batch, transfer_cond);
}

template <typename T, typename K>
void QueueUtil::push(SyclQueue& Q, FullereneBatch<T, K>& dst_batch, FullereneQueue<T, K>& src_queue, ConditionFunctor transfer_cond){
    push_impl(Q, dst_batch, src_queue, transfer_cond);
}

template <typename DstBatch, typename SrcBatch>
void QueueUtil::push_impl(SyclQueue& Q, DstBatch& dst_batch, SrcBatch& src_batch, ConditionFunctor condition){
if(dst_batch.N_ != src_batch.N_ || src_batch.Nf_ != dst_batch.Nf_) throw std::runtime_error("Fullerenes have different sizes");
    constexpr bool dst_is_queue = is_fullerene_queue<DstBatch>::value;
    //If we are pushing from a FullereneQueue, we need to compute the valid indices of the destination batch
    //The queue is assumed to be contiguous, so we know that all fullernes between the front and back are valid
    //Vice versa, if we are pushing to a FullereneQueue, we need to compute the valid indices of the source batch
    auto& indices = dst_is_queue ? src_batch.m_.valid_indices_ : dst_batch.m_.valid_indices_;
    auto& flags = dst_is_queue ? src_batch.m_.flags_ : dst_batch.m_.flags_;
    using int_t = std::decay_t<decltype(*indices.begin())>;
    using float_t = std::decay_t<decltype(src_batch.d_.X_cubic_[0][0])>;
    int_t init = 0;
    std::transform_exclusive_scan(std::execution::par_unseq,
        flags.begin(), flags.end(),
        indices.begin(),
        init,
        Plus{},
        [condition](auto f){return static_cast<int_t>(condition(f));}
    );
    Q.wait();
    
    auto src = FullereneBatchView(src_batch, 0, src_batch.capacity());
    auto dst = FullereneBatchView(dst_batch, 0, dst_batch.capacity());

    auto num_valid = condition(flags.back()) + indices.back(); // Number of valid fullerenes in the source batch
    if constexpr (dst_is_queue) dst_batch.resize( num_valid + dst_batch.size()  );
    auto num_transfer = dst_is_queue ? num_valid : std::min((int_t)num_valid, (int_t)src_batch.size());

    FullereneQueue<float_t, int_t>* queue;
    if constexpr (dst_is_queue) {queue = &dst_batch;} else {queue = &src_batch;}


    auto& queue_front = queue->front_;
    auto& queue_back = queue->back_;
    auto& queue_size = queue->size_;
    auto queue_capacity = queue->capacity();

    if (dst_is_queue && queue_front < 0) {queue_front = 0; queue_back = 0;}

    std::for_each(indices.begin(), indices.end(), [&](auto& it){
        auto index_in_indices = &it - &(*indices.begin());
        if (condition(flags[index_in_indices]) && (dst_is_queue || index_in_indices < num_transfer)){ 
            auto index_in_dst = (dst_is_queue ? (queue_back + it + 1) : index_in_indices) % dst_batch.capacity();
            auto index_in_src = (dst_is_queue ? index_in_indices : queue_front + it) % src_batch.capacity();
            decltype(dst[0])::copy(dst[index_in_dst], src[index_in_src]);
        }
    });

    /* primitives::for_each(Q, indices, [=](auto& it){
        auto index_in_indices = &it - &(*indices.begin());
        if (condition(flags[index_in_indices]) && (dst_is_queue || index_in_indices < num_transfer)){ 
            auto index_in_dst = (dst_is_queue ? (queue_back + it + 1) : index_in_indices) % dst_batch.capacity();
            auto index_in_src = (dst_is_queue ? index_in_indices : queue_front + it) % src_batch.capacity();
            decltype(dst[0])::copy(dst[index_in_dst], src[index_in_src]);
        }
    }); */
    

    queue_size += dst_is_queue ? (num_transfer) : (-num_transfer);
    if constexpr(dst_is_queue) { queue_back = (queue_back + num_transfer) % queue_capacity;} else { queue_front = (queue_front + num_transfer) % queue_capacity;}
    if (queue_size == 0) {queue_front = -1; queue_back = -1;}
    if (!dst_is_queue) dst_batch.size_ = std::min((int)dst_batch.capacity(), (int)(dst_batch.size() + num_transfer));
}

template <typename T, typename K>
Fullerene<T,K> FullereneBatch<T,K>::operator[](size_t index) const{
    if(index >= capacity_) {
        std::cout << "Index: " << index << " Size: " << capacity_ << std::endl;
        throw OOR_ERROR(index, capacity_);}
    return FullereneBatchView(*this)[index];
}

template struct SyclVector<int8_t>;
template struct SyclVector<int16_t>;
template struct SyclVector<int32_t>;
template struct SyclVector<int64_t>;
template struct SyclVector<uint8_t>;
template struct SyclVector<uint16_t>;
template struct SyclVector<uint32_t>;
template struct SyclVector<uint64_t>;
template struct SyclVector<float>;
template struct SyclVector<double>;
template struct SyclVector<std::byte>;
template struct SyclVector<StatusFlag>;
template struct SyclVector<bool>;
template struct SyclVector<std::array<double,2>>;
template struct SyclVector<std::array<double,3>>;
template struct SyclVector<std::array<float,2>>;
template struct SyclVector<std::array<float,3>>;
template struct SyclVector<std::array<uint16_t,3>>;
template struct SyclVector<std::array<uint32_t,3>>;
template struct SyclVector<std::array<uint16_t,6>>;
template struct SyclVector<std::array<uint32_t,6>>;
//template struct SyclVector<NodeNeighbours<uint16_t>>;
//template struct SyclVector<NodeNeighbours<uint32_t>>;
//template struct SyclVector<Constants<float,uint16_t>>;

template struct Span<int8_t>;
template struct Span<int16_t>;
template struct Span<int32_t>;
template struct Span<int64_t>;
template struct Span<uint8_t>;
template struct Span<uint16_t>;
template struct Span<uint32_t>;
template struct Span<uint64_t>;
template struct Span<float>;
template struct Span<double>;
template struct Span<std::byte>;
template struct Span<StatusFlag>;
template struct Span<bool>;
template struct Span<std::array<double,2>>;
template struct Span<std::array<double,3>>;
template struct Span<std::array<float,2>>;
template struct Span<std::array<float,3>>;
template struct Span<std::array<uint16_t,3>>;
template struct Span<std::array<uint32_t,3>>;
template struct Span<std::array<uint16_t,6>>;
template struct Span<std::array<uint32_t,6>>;
//template struct Span<NodeNeighbours<uint16_t>>;
//template struct Span<NodeNeighbours<uint32_t>>;
//template struct Span<Constants<float,uint16_t>>;


template struct Fullerene<float,uint16_t>;
template struct Fullerene<double,uint16_t>;
template struct Fullerene<float,uint32_t>;
template struct Fullerene<double,uint32_t>;

template struct FullereneBatch<float,uint16_t>;
template struct FullereneBatch<double,uint16_t>;
template struct FullereneBatch<float,uint32_t>;
template struct FullereneBatch<double,uint32_t>;

template struct FullereneQueue<float,uint16_t>;
template struct FullereneQueue<double,uint16_t>;
template struct FullereneQueue<float,uint32_t>;
template struct FullereneQueue<double,uint32_t>;

template std::ostream& operator<<(std::ostream& os, const SyclVector<float>& vec);
template std::ostream& operator<<(std::ostream& os, const SyclVector<double>& vec);
template std::ostream& operator<<(std::ostream& os, const SyclVector<uint16_t>& vec);
template std::ostream& operator<<(std::ostream& os, const SyclVector<uint32_t>& vec);
template std::ostream& operator<<(std::ostream& os, const SyclVector<int>& vec);
template std::ostream& operator<<(std::ostream& os, const SyclVector<size_t>& vec);
template std::ostream& operator<<(std::ostream& os, const SyclVector<StatusFlag>& vec);
template std::ostream& operator<<(std::ostream& os, const SyclVector<bool>& vec);
template std::ostream& operator<<(std::ostream& os, const SyclVector<std::array<double,2>>& vec);
template std::ostream& operator<<(std::ostream& os, const SyclVector<std::array<double,3>>& vec);
template std::ostream& operator<<(std::ostream& os, const SyclVector<std::array<float,2>>& vec);
template std::ostream& operator<<(std::ostream& os, const SyclVector<std::array<float,3>>& vec);
template std::ostream& operator<<(std::ostream& os, const SyclVector<std::array<uint16_t,3>>& vec);
template std::ostream& operator<<(std::ostream& os, const SyclVector<std::array<uint32_t,3>>& vec);
template std::ostream& operator<<(std::ostream& os, const SyclVector<std::array<uint16_t,6>>& vec);
template std::ostream& operator<<(std::ostream& os, const SyclVector<std::array<uint32_t,6>>& vec);


template std::ostream& operator<<(std::ostream& os, const Span<float>& vec);
template std::ostream& operator<<(std::ostream& os, const Span<double>& vec);
template std::ostream& operator<<(std::ostream& os, const Span<uint16_t>& vec);
template std::ostream& operator<<(std::ostream& os, const Span<uint32_t>& vec);

template std::ostream& operator<<(std::ostream& os, const Span<int>& vec);
template std::ostream& operator<<(std::ostream& os, const Span<size_t>& vec);
template std::ostream& operator<<(std::ostream& os, const Span<StatusFlag>& vec);
template std::ostream& operator<<(std::ostream& os, const Span<bool>& vec);
template std::ostream& operator<<(std::ostream& os, const Span<std::array<double,2>>& vec);
template std::ostream& operator<<(std::ostream& os, const Span<std::array<double,3>>& vec);
template std::ostream& operator<<(std::ostream& os, const Span<std::array<float,2>>& vec);
template std::ostream& operator<<(std::ostream& os, const Span<std::array<float,3>>& vec);
template std::ostream& operator<<(std::ostream& os, const Span<std::array<uint16_t,3>>& vec);
template std::ostream& operator<<(std::ostream& os, const Span<std::array<uint32_t,3>>& vec);
template std::ostream& operator<<(std::ostream& os, const Span<std::array<uint16_t,6>>& vec);
template std::ostream& operator<<(std::ostream& os, const Span<std::array<uint32_t,6>>& vec);

template std::ostream& operator<<(std::ostream& os, const Fullerene<float,uint16_t>& vec);
template std::ostream& operator<<(std::ostream& os, const Fullerene<double,uint16_t>& vec);
template std::ostream& operator<<(std::ostream& os, const Fullerene<float,uint32_t>& vec);
template std::ostream& operator<<(std::ostream& os, const Fullerene<double,uint32_t>& vec);

template std::ostream& operator<<(std::ostream& os, const FullereneBatch<float,uint16_t>& vec);
template std::ostream& operator<<(std::ostream& os, const FullereneBatch<double,uint16_t>& vec);
template std::ostream& operator<<(std::ostream& os, const FullereneBatch<float,uint32_t>& vec);
template std::ostream& operator<<(std::ostream& os, const FullereneBatch<double,uint32_t>& vec);

template std::ostream& operator<<(std::ostream& os, const FullereneBatchView<float,uint16_t>& vec);
template std::ostream& operator<<(std::ostream& os, const FullereneBatchView<double,uint16_t>& vec);
template std::ostream& operator<<(std::ostream& os, const FullereneBatchView<float,uint32_t>& vec);
template std::ostream& operator<<(std::ostream& os, const FullereneBatchView<double,uint32_t>& vec);

template std::ostream& operator<<(std::ostream& os, const ReferenceWrapper<int>& ref);
template std::ostream& operator<<(std::ostream& os, const ReferenceWrapper<size_t>& ref);
template std::ostream& operator<<(std::ostream& os, const ReferenceWrapper<float>& ref);
template std::ostream& operator<<(std::ostream& os, const ReferenceWrapper<double>& ref);
template std::ostream& operator<<(std::ostream& os, const ReferenceWrapper<uint16_t>& ref);
template std::ostream& operator<<(std::ostream& os, const ReferenceWrapper<StatusFlag>& ref);



template void QueueUtil::push(SyclQueue& Q, FullereneQueue<float,uint16_t>& dst_queue, FullereneBatch<float,uint16_t>& src_batch, ConditionFunctor transfer_cond);
template void QueueUtil::push(SyclQueue& Q, FullereneBatch<float,uint16_t>& dst_batch, FullereneQueue<float,uint16_t>& src_queue, ConditionFunctor transfer_cond);
template void QueueUtil::push(SyclQueue& Q, FullereneQueue<float,uint32_t>& dst_queue, FullereneBatch<float,uint32_t>& src_batch, ConditionFunctor transfer_cond);
template void QueueUtil::push(SyclQueue& Q, FullereneBatch<float,uint32_t>& dst_batch, FullereneQueue<float,uint32_t>& src_queue, ConditionFunctor transfer_cond);
//template void push_impl (SyclQueue& Q, FullereneQueue<float,uint16_t>& dst_queue, FullereneBatch<float,uint16_t>& src_batch, ConditionFunctor condition);
