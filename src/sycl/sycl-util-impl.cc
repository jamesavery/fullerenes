#include <fullerenes/sycl-wrappers.hh>
#include <CL/sycl.hpp>
#include <unordered_map>
#include <exception>
#include <limits>
#include "primitives.hh"
#ifndef DEVICE_CAST
    #define DEVICE_CAST(x,ix) (reinterpret_cast<const sycl::device*>(x)[ix])
#endif
#include "queue-impl.cc"
//#define HOST_ALLOCATOR_QUEUE reinterpret_cast<sycl::queue*>(QueueWrapper::allocator_queue_.data())[0]

using namespace sycl;

std::vector<std::byte> queue_to_bytes(const sycl::queue& q){
    return std::vector<std::byte>(reinterpret_cast<const std::byte*>(&q), reinterpret_cast<const std::byte*>(&q) + sizeof(sycl::queue));
}

std::vector<std::byte> device_vector_to_bytes(const std::vector<sycl::device>& devices){
    return std::vector<std::byte>(reinterpret_cast<const std::byte*>(devices.data()), reinterpret_cast<const std::byte*>(devices.data()) + sizeof(sycl::device)*devices.size());
}

template <typename T>
SyclVector<T>::SyclVector() : size_(0), capacity_(0), data_(nullptr) {}

template <typename T>
SyclVector<T>::SyclVector(size_t size) : size_(size), capacity_(size) {
    static sycl::queue Q(default_selector_v, sycl::property::queue::in_order{});
    std::cout << "Allocating " << size << " elements" << std::endl;
    data_ = sycl::malloc_shared<T>(size, Q);
    Q.wait();
}

template <typename T>
SyclVector<T>::SyclVector(size_t size, T value) : size_(size), capacity_(size) {
    static sycl::queue Q(default_selector_v, sycl::property::queue::in_order{});
    data_ = sycl::malloc_shared<T>(size, Q);
    for(size_t i = 0; i < size; i++) data_[i] = value;
    Q.wait();
}

template <typename T>
SyclVector<T>::SyclVector(const SyclVector<T>& other) : size_(other.size_), capacity_(other.capacity_) {
    static sycl::queue Q(default_selector_v, sycl::property::queue::in_order{});
    data_ = sycl::malloc_shared<T>(capacity_, Q);
    for(size_t i = 0; i < size_; i++) data_[i] = other.data_[i];
    Q.wait();
}
template <typename T>
SyclVector<T>::SyclVector(SyclVector<T>&& other) : size_(other.size_), capacity_(other.capacity_), data_(other.data_) {
    other.data_ = nullptr;
}

template <typename T>
SyclVector<T>::~SyclVector() {static sycl::queue Q(default_selector_v, sycl::property::queue::in_order{}); if(data_) sycl::free(data_, Q);}

template <typename T>
void SyclVector<T>::fill(T data) {
    for(size_t i = 0; i < size_; i++) data_[i] = data;
}

template <typename T>
T*  SyclVector<T>::data() const {return data_;}

template <typename T>
size_t SyclVector<T>::size() const {return size_;}

template <typename T>
size_t SyclVector<T>::capacity() const {return capacity_;}

template <typename T>
void SyclVector<T>::resize(size_t new_size) {
    static sycl::queue Q(default_selector_v, sycl::property::queue::in_order{});
    if(new_size > capacity_){
        T* new_data = sycl::malloc_shared<T>(new_size, Q);
        memcpy(new_data, data_, size_*sizeof(T));
        Q.wait();
        if(data_) sycl::free(data_, Q);
        data_ = new_data;
        capacity_ = new_size;
    }
    size_ = new_size;
}

template <typename T>
void SyclVector<T>::reserve(size_t new_capacity) {
    static sycl::queue Q(default_selector_v, sycl::property::queue::in_order{});
    if(new_capacity > capacity_){
        T* new_data = sycl::malloc_shared<T>(new_capacity, Q);
        memcpy(new_data, data_, size_*sizeof(T));
        if(data_) sycl::free(data_, Q);
        data_ = new_data;
        capacity_ = new_capacity;
    }
}

template <typename T>
void SyclVector<T>::clear() {size_ = 0;}

template <typename T>
T& SyclVector<T>::operator[](size_t index) {return data_[index];}

template <typename T>
const T& SyclVector<T>::operator[](size_t index) const {return data_[index];}

template <typename T>
T& SyclVector<T>::at(size_t index) {
    if(index >= size_) throw std::out_of_range("Index out of range");
    return data_[index];
}

template <typename T>
const T& SyclVector<T>::at(size_t index) const {
    if(index >= size_) throw std::out_of_range("Index out of range");
    return data_[index];
}

template <typename T>
void SyclVector<T>::push_back(const T& value) {
    if(size_ == capacity_){
        size_t new_capacity = capacity_ == 0 ? 1 : 2*capacity_;
        reserve(new_capacity);
    }
    data_[size_++] = value;
}

template <typename T>
void SyclVector<T>::pop_back() {
    if(size_ > 0) size_--;
}

template <typename T>
SyclVector<T>& SyclVector<T>::operator=(const SyclVector<T>& other) {
    static sycl::queue Q(default_selector_v, sycl::property::queue::in_order{});
    if(data_) sycl::free(data_, Q);
    size_ = other.size_;
    capacity_ = other.capacity_;
    data_ = sycl::malloc_shared<T>(capacity_, Q);
    for(size_t i = 0; i < size_; i++) data_[i] = other.data_[i];
    Q.wait();
    return *this;
}

template <typename T>
SyclVector<T>& SyclVector<T>::operator=(SyclVector<T>&& other){
    data_ = other.data_;
    size_ = other.size_;
    capacity_ = other.capacity_;
    other.data_ = nullptr;
}

template <typename T>
bool SyclVector<T>::operator==(const SyclVector<T>& other) const {
    if(size_ != other.size_) return false;
    if constexpr (is_floating_point_v<T>){
        float eps = std::numeric_limits<T>::epsilon() * 10;
        for(size_t i = 0; i < size_; i++){
            T max_v = std::max<T>(std::abs<T>(data_[i]), std::abs<T>(other.data_[i]));
            if(std::abs<T>(data_[i] - other.data_[i]) / (max_v > eps ? max_v : 1) > eps) return false;
        }
    } else{
        return std::equal(begin(), end(), other.begin());
    }
}

template <typename U>
std::ostream& operator<<(std::ostream& os, const SyclVector<U>& vec) {
    os << "[";
    for(size_t i = 0; i < vec.size_; i++){
        os << vec.data_[i];
        if(i < vec.size_ - 1) os << ", ";
    }
    os << "]\n";
    return os;
}



template <typename T>
T* SyclVector<T>::begin() const {return (T*) (data_);}

template <typename T>
T* SyclVector<T>::end() const {return (T*) (data_ + size_);}

template <typename T>
SyclVector<T>::operator Span<T>() const {return Span<T>(data_, size_);}

template <typename T>
Span<T> SyclVector<T>::operator()(size_t offset, size_t count) const {
    assert(offset + count <= size_ && "Span out of range of vector");
    return Span<T>(data_ + offset, count);
}

template <typename T>
bool Span<T>::operator==(const Span<T>& other) const {
    if(size_ != other.size_) return false;
    if(data() == other.data()) return true;
    if constexpr (is_floating_point_v<T>){
        float eps = std::numeric_limits<T>::epsilon() * 10;
        for(size_t i = 0; i < size_; i++){
            T max_v = std::max<T>(std::abs<T>(data_[i]), std::abs<T>(other.data_[i]));
            if(std::abs<T>(data_[i] - other.data_[i]) / (max_v > eps ? max_v : 1) > eps) return false;
        }
    } else{
        return std::equal(begin(), end(), other.begin());
    }
}

template <typename T, typename K>
Fullerene<T,K>::Fullerene(const FullereneDataMembers<Span, T,K>& data, size_t N, size_t Nf) : N_(N), Nf_(Nf) {
    d_ = data;
}


template <typename T, typename K>
bool Fullerene<T,K>::operator==(const Fullerene<T,K>& other) const {
    if(N_ != other.N_ || Nf_ != other.Nf_) return false;
    return std::apply([&](auto&&... args){return ((std::get<0>(args) == std::get<1>(args)) && ...);}, forward_merge_tuples(d_.to_tuple(),other.d_.to_tuple()));
}

template <typename T, typename K>
FullereneBatch<T,K>::FullereneBatch() : N_(0), Nf_(0), size_(0), capacity_(0), front_(-1), back_(-1) {}


template <size_t N, typename... Args, std::size_t... I>
void resize_all(std::array<int,N>&& sizes, std::tuple<Args...>&& args, std::index_sequence<I...>){
    (std::get<I>(args).resize(sizes[I]), ...);
}

template <typename T, typename K>
FullereneBatch<T,K>::FullereneBatch(size_t N, int capacity) :
    N_(N),
    Nf_(N/2 + 2),
    capacity_(capacity),
    size_(0),
    front_(-1),
    back_(-1)
{
    resize_all(d_.get_size_factors(N_, capacity_), d_.to_tuple(), std::make_index_sequence<std::tuple_size_v<decltype(d_.to_tuple())>>{});
    resize_all(m_.get_size_factors(capacity_), m_.to_tuple(), std::make_index_sequence<std::tuple_size_v<decltype(m_.to_tuple())>>{});
}

template <typename T, typename K>
void FullereneBatch<T,K>::resize(size_t new_capacity) {
    if(new_capacity == capacity_) return;
    capacity_ = new_capacity;
    resize_all(d_.get_size_factors(N_, capacity_), d_.to_tuple(), std::make_index_sequence<std::tuple_size_v<decltype(d_.to_tuple())>>{});
    resize_all(m_.get_size_factors(capacity_), m_.to_tuple(), std::make_index_sequence<std::tuple_size_v<decltype(m_.to_tuple())>>{});
}

template <typename T, typename K>
void FullereneBatch<T,K>::prepare_for_push_back(const neighbours_t& neighbours, bool is_cubic){
    if((is_cubic && N_ != neighbours.size()) || (!is_cubic && Nf_ != neighbours.size())) throw std::runtime_error("Graph size is incompatible with the batch");
    if(capacity_ == size_) resize(capacity_ > 0 ? 2*capacity_ : 1);
    if(back_ == -1) {back_ = 0; front_ = 0;}
    else {back_ = (back_ + 1) % capacity_;}
    size_++;
}

template <typename T, typename K>
void FullereneBatch<T,K>::push_back(const neighbours_t& neighbours, bool is_cubic, const int ID) {
    prepare_for_push_back(neighbours, is_cubic);
    size_t N = is_cubic ? N_ : Nf_;
    for(size_t i = 0; i < N; i++){
        if(!is_cubic) d_.deg_[i + back_*Nf_] = neighbours[i].size();
        for(size_t j = 0; j < neighbours[i].size(); j++){
            if(is_cubic){d_.A_cubic_[i*3 + j + back_*N_*3] = neighbours[i][j];}
            else{d_.A_dual_[i*6 + j + back_*Nf_*6] = neighbours[i][j];}
        }
    }
    m_.ID_[back_] = ID;
    m_.flags_[back_] = is_cubic ? StatusFlag::CUBIC_INITIALIZED : StatusFlag::DUAL_INITIALIZED;
    m_.iterations_[back_] = 0;
}

template <typename T, typename K>
void FullereneBatch<T,K>::push_back(const Graph& G, const int ID) {push_back(G.neighbours, false, ID);}

template <typename T, typename K>
void FullereneBatch<T,K>::push_back(const Polyhedron& P, const int ID) {
    push_back(P.neighbours, true, ID);
    //If floating point types are not the same, we need to convert
    if(!std::is_same<T, decltype(P.points[0][0])>::value){
        
    for(size_t i = 0; i < N_; i++){
        for(size_t j = 0; j < 3; j++){
            d_.X_cubic_[i*3 + j + back_*N_*3] = P.points[i][j];
        }
    } //Otherwise a contiguous memory copy will do
    } else{ memcpy(d_.X_cubic_.data() + back_*N_*3, P.points.data(), N_*3*sizeof(T));}
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
template <typename T, typename K> 
void push_impl(SyclContext& ctx, size_t Qid, FullereneBatch<T,K>& dst_batch, FullereneBatch<T,K>& src_batch, bool dst_is_queue, ConditionFunctor condition) {
    if(dst_batch.N_ != src_batch.N_ || dst_batch.Nf_ != dst_batch.Nf_) throw std::runtime_error("Fullerenes have different sizes");
    auto& Q = reinterpret_cast<sycl::queue*>(ctx.queues_.data())[Qid];

    //If we are pushing from a FullereneQueue, we need to compute the valid indices of the destination batch
    //The queue is assumed to be contiguous, so we know that all fullernes between the front and back are valid
    auto& indices = dst_is_queue ? src_batch.m_.valid_indices_ : dst_batch.m_.valid_indices_;
    auto& flags = dst_is_queue ? src_batch.m_.flags_ : dst_batch.m_.flags_;
    primitives::transform_exclusive_scan(Q, 
        flags.begin(), flags.end(),
        indices.begin(),
        K(0),
        std::plus<K>(),
        condition
    );

    primitives::transform(Q, 
        flags.begin(), flags.end(),
        indices.begin(),
        [=](K x){return x == 0 ? 0 : 1;}
    );

    //Trivially copyable views of the batches (Necessary property for SYCL kernels)
    FullereneBatchView<T,K> src(src_batch); 
    FullereneBatchView<T,K> dst(dst_batch);

    //The exclusive scan will contain the number of valid indices at the back
    size_t num_valid = indices.back();
    int diff = num_valid - (dst_batch.capacity() - dst_batch.size());
    if(dst_is_queue && diff > 0) dst_batch.resize(dst_batch.capacity() + diff);
    size_t num_transfer = dst_is_queue ? num_valid : src_batch.size();
    int queue_front = dst_is_queue ? dst_batch.front() : src_batch.front();


    Q.submit([&](handler& cgh){
        cgh.parallel_for<FullereneBatchFill<T,K>>(nd_range<1>(range{}, range{1}), [=](nd_item<1> item){
            size_t bid = item.get_group_linear_id();
            if(!condition(dst_is_queue ? src.m_.flags_[bid] : dst.m_.flags_[bid])) return;

            //If we are pushing from a queue, the destination index is determined by the valid indices
            //If we are pushing onto a queue, the destination index is determined by the queue front
            //Similarly logic applies to the source index
            auto i = dst_is_queue ? queue_front + bid : dst.m_.valid_indices_[bid]; 
            auto j = dst_is_queue ? src.m_.valid_indices_[bid] : queue_front + bid;



            auto grp = item.get_group();
            auto dst_fullerene = dst[i];
            auto src_fullerene = src[j];
            auto span_pairs = forward_merge_tuples(src_fullerene.d_.to_tuple(), dst_fullerene.d_.to_tuple());
            std::apply([&grp](auto&&... args){copy_spans(grp, args...);}, span_pairs);
            dst.m_.flags_[i] = src.m_.flags_[j];
            dst.m_.ID_[i] = src.m_.ID_[j];
            dst.m_.iterations_[i] = src.m_.iterations_[j];
        });
    });
    dst_batch.size_ += num_transfer;
    src_batch.size_ -= num_transfer;
}

/* template <typename T, typename K>
void FullereneQueue<T,K>::push(SyclContext& ctx, size_t Qid, FullereneBatch<T,K>& batch, ConditionFunctor condition) {
    push_impl(ctx, Qid, *this, batch, condition, condition);
}
 */
template <typename T, typename K>
Fullerene<T,K> FullereneBatch<T,K>::operator[](size_t index) const{
    if(index >= size_) throw std::out_of_range("Index out of range");
    auto own_data = d_.to_tuple();
    auto size_data = d_.get_size_factors(N_, capacity_);
    auto seq = std::make_index_sequence<std::tuple_size_v<decltype(d_.to_tuple())>>{};
    
    auto span_tuple = fullerene_detail::with_sequence<std::tuple_size_v<decltype(d_.to_tuple())>>(
        [&](auto... I){
            return std::make_tuple(Span(std::get<I>(own_data).data() + size_data[I] * index, size_data[I])...);
        }
    );
    FullereneDataMembers<Span, T,K> data;
    auto fdata_tuple = data.to_tuple();

    fullerene_detail::construct_spans(index, std::move(size_data), fdata_tuple, own_data);
    
    return Fullerene<T,K>(data,N_,Nf_);
}

SyclContext::SyclContext() : 
    cpus_(device_vector_to_bytes(sycl::device::get_devices(sycl::info::device_type::cpu))), 
    gpus_(device_vector_to_bytes(sycl::device::get_devices(sycl::info::device_type::gpu))),
    accelerators_(device_vector_to_bytes(sycl::device::get_devices(sycl::info::device_type::accelerator)))
    {}

SyclContext::~SyclContext() = default;

std::string SyclContext::device_get_name(Device dev) const  {
    switch (dev.second){
        case DeviceType::CPU:
            return DEVICE_CAST(cpus_.data(), dev.first).get_info<sycl::info::device::name>();
        case DeviceType::GPU:
            return DEVICE_CAST(gpus_.data(), dev.first).get_info<sycl::info::device::name>();
        case DeviceType::ACCELERATOR:
            return DEVICE_CAST(accelerators_.data(), dev.first).get_info<sycl::info::device::name>();
        default:
            return "Unknown Device";
    }
}

size_t SyclContext::device_get_count(DeviceType type) const{
    switch(int(type)){
        case 0: return cpus_.size()/sizeof(sycl::device);
        case 1: return gpus_.size()/sizeof(sycl::device);
        case 2: return accelerators_.size()/sizeof(sycl::device);
        default: return 0;
    }
}

size_t SyclContext::device_get_count() const {
    return (cpus_.size() + gpus_.size() + accelerators_.size())/sizeof(sycl::device);
}
    
size_t SyclContext::device_get_property(Device dev, DeviceProperty property) const {
    if(dev.first >= device_get_count(dev.second)) throw std::out_of_range("Device index out of range");
    const sycl::device& d = dev.second == DeviceType::CPU ? DEVICE_CAST(cpus_.data(), dev.first) : dev.second == DeviceType::GPU ? DEVICE_CAST(gpus_.data(), dev.first) : DEVICE_CAST(accelerators_.data(), dev.first);

    switch(int(property)){
        case 0: return d.get_info<sycl::info::device::max_work_group_size>();
        case 1: return d.get_info<sycl::info::device::max_clock_frequency>();
        case 2: return d.get_info<sycl::info::device::max_compute_units>();
        case 3: return d.get_info<sycl::info::device::max_mem_alloc_size>();
        case 4: return d.get_info<sycl::info::device::global_mem_size>();
        case 5: return d.get_info<sycl::info::device::local_mem_size>();
        case 7: return d.get_info<sycl::info::device::max_constant_args>();
        default: std::cerr << "Unknown property" << std::endl; return 0;
    }
}

void SyclContext::queue_push_back(Device dev, bool in_order) {
    if(dev.first >= device_get_count(dev.second)) throw std::out_of_range("Device index out of range");
    queues_.resize(queues_.size() + sizeof(sycl::queue));
    const void* device_arrays[] = { cpus_.data(), gpus_.data(), accelerators_.data() };
    reinterpret_cast<sycl::queue*>(queues_.data())[queues_.size() - 1] = in_order 
        ? sycl::queue(DEVICE_CAST(device_arrays[(int)dev.second], dev.first), sycl::property::queue::in_order{}) 
        : sycl::queue(DEVICE_CAST(device_arrays[(int)dev.second], dev.first));
}

SyclQueue::SyclQueue() : device_({0, DeviceType::CPU}), in_order_(true) {
    impl_ = std::make_unique<SyclQueueImpl>(device_, in_order_);
}

SyclQueue::SyclQueue(Device dev, bool in_order) : device_(dev), in_order_(in_order) {
    impl_ = std::make_unique<SyclQueueImpl>(dev, in_order);
}

SyclQueue::~SyclQueue() = default;

void SyclQueue::wait() const {impl_->wait();}
void SyclQueue::wait_and_throw() const {impl_->wait_and_throw();}


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

template struct Fullerene<float,uint16_t>;
template struct Fullerene<double,uint16_t>;

template struct FullereneBatch<float,uint16_t>;
template struct FullereneBatch<double,uint16_t>;

template std::ostream& operator<<(std::ostream& os, const SyclVector<float>& vec);
template std::ostream& operator<<(std::ostream& os, const SyclVector<double>& vec);
template std::ostream& operator<<(std::ostream& os, const SyclVector<uint16_t>& vec);
template std::ostream& operator<<(std::ostream& os, const SyclVector<int>& vec);
template std::ostream& operator<<(std::ostream& os, const SyclVector<size_t>& vec);