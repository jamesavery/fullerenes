#include <fullerenes/sycl-wrappers.hh>
#include <CL/sycl.hpp>
#include <unordered_map>
#define DEVICE_CAST(x) (*static_cast<const sycl::device*>(x))
//#define HOST_ALLOCATOR_QUEUE reinterpret_cast<sycl::queue*>(QueueWrapper::allocator_queue_.data())[0]

using namespace sycl;


std::vector<std::byte> queue_to_bytes(const sycl::queue& q){
    return std::vector<std::byte>(reinterpret_cast<const std::byte*>(&q), reinterpret_cast<const std::byte*>(&q) + sizeof(sycl::queue));
}

std::vector<std::byte> device_vector_to_bytes(const std::vector<sycl::device>& devices){
    return std::vector<std::byte>(reinterpret_cast<const std::byte*>(devices.data()), reinterpret_cast<const std::byte*>(devices.data()) + sizeof(sycl::device)*devices.size());
}

std::vector<std::byte> DeviceWrapper::__FULLERENE_ALLACCELERATORS__  = device_vector_to_bytes(device::get_devices(info::device_type::accelerator));
std::vector<std::byte> DeviceWrapper::__FULLERENE_ALLGPUDEVICES__ = device_vector_to_bytes(device::get_devices(info::device_type::gpu));
std::vector<std::byte> DeviceWrapper::__FULLERENE_ALLCPUDEVICES__ = device_vector_to_bytes(device::get_devices(info::device_type::cpu));

//std::vector<std::byte> QueueWrapper::allocator_queue_ = queue_to_bytes(sycl::queue(sycl::cpu_selector{}));

static sycl::queue HOST_ALLOCATOR_QUEUE = sycl::queue(sycl::cpu_selector{});

template <typename T>
SyclVector<T>::SyclVector() : size_(0), capacity_(0), data_(nullptr) {}

template <typename T>
SyclVector<T>::SyclVector(size_t size) : size_(size), capacity_(size), data_(nullptr) {
    std::cout << "Allocating " << size << " elements" << std::endl;
    data_ = sycl::malloc_shared<T>(size, HOST_ALLOCATOR_QUEUE);
    HOST_ALLOCATOR_QUEUE.wait();
}

template <typename T>
SyclVector<T>::SyclVector(size_t size, T value) : size_(size), capacity_(size) {
    data_ = sycl::malloc_shared<T>(size, HOST_ALLOCATOR_QUEUE);
    for(size_t i = 0; i < size; i++) data_[i] = value;
    HOST_ALLOCATOR_QUEUE.wait();
}

template <typename T>
SyclVector<T>::SyclVector(const SyclVector<T>& other) : size_(other.size_), capacity_(other.capacity_) {
    data_ = sycl::malloc_shared<T>(capacity_, HOST_ALLOCATOR_QUEUE);
    for(size_t i = 0; i < size_; i++) data_[i] = other.data_[i];
    HOST_ALLOCATOR_QUEUE.wait();
}
template <typename T>
SyclVector<T>::SyclVector(SyclVector<T>&& other) : size_(other.size_), capacity_(other.capacity_), data_(other.data_) {
    other.data_ = nullptr;
}

template <typename T>
SyclVector<T>::~SyclVector() {if(data_) sycl::free(data_, HOST_ALLOCATOR_QUEUE);}

template <typename T>
void SyclVector<T>::fill(T data) {
    for(size_t i = 0; i < size_; i++) data_[i] = data;
}

template <typename T>
T*  SyclVector<T>::data() {return data_;}

template <typename T>
size_t SyclVector<T>::size() {return size_;}

template <typename T>
size_t SyclVector<T>::capacity() {return capacity_;}

template <typename T>
void SyclVector<T>::resize(size_t new_size) {
    if(new_size > capacity_){
        T* new_data = sycl::malloc_shared<T>(new_size, HOST_ALLOCATOR_QUEUE);
        for(size_t i = 0; i < size_; i++) new_data[i] = data_[i];
        HOST_ALLOCATOR_QUEUE.wait();
        if(data_) sycl::free(data_, HOST_ALLOCATOR_QUEUE);
        data_ = new_data;
        capacity_ = new_size;
    }
    size_ = new_size;
}

template <typename T>
void SyclVector<T>::reserve(size_t new_capacity) {
    if(new_capacity > capacity_){
        T* new_data = sycl::malloc_shared<T>(new_capacity, HOST_ALLOCATOR_QUEUE);
        for(size_t i = 0; i < size_; i++) new_data[i] = data_[i];
        if(data_) sycl::free(data_, HOST_ALLOCATOR_QUEUE);
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
    if(data_) sycl::free(data_, HOST_ALLOCATOR_QUEUE);
    size_ = other.size_;
    capacity_ = other.capacity_;
    data_ = sycl::malloc_shared<T>(capacity_, HOST_ALLOCATOR_QUEUE);
    for(size_t i = 0; i < size_; i++) data_[i] = other.data_[i];
    HOST_ALLOCATOR_QUEUE.wait();
    return *this;
}


template <typename T, typename K>
FullereneIsomer<T,K>::FullereneIsomer() : N_(0), Nf_(0) {}

template <typename T, typename K>
FullereneIsomer<T,K>::~FullereneIsomer() = default;

template <typename T, typename K>
FullereneIsomer<T,K>::FullereneIsomer(size_t N, bool allocate_hessian_etc) : 
    N_(N),
    Nf_(N/2 + 2),
    faces_(Nf_*6), 
    deg_(Nf_),
    X_cubic_(N_*3),
    X_dual_(Nf_*6),
    A_cubic_(N_*3),
    A_dual_(Nf_*6),
    gradient_(N_*3)
{

    if(allocate_hessian_etc){
        hessian_ = SyclVector<T>(N_*90);
        eigenvalues_ = SyclVector<T>(N_*3);
        eigenvectors_ = SyclVector<T>(N_*3*N_*3);
    }
}

template <typename T, typename K>
FullereneIsomer<T,K>::FullereneIsomer(const neighbours_t& G, bool is_cubic, bool allocate_hessian_etc)
{   
    N_ = (is_cubic ? G.size() : (G.size()-2)*2);
    Nf_ = (N_/2 + 2);
    faces_.resize(Nf_*6);
    deg_.resize(Nf_);
    X_cubic_.resize(N_*3);
    X_dual_.resize(Nf_*6);
    A_cubic_.resize(N_*3);
    A_dual_.resize(Nf_*6);
    gradient_.resize(N_*3);
    assert(N_ > 0 && Nf_ > 0 && deg_.capacity() > 0);
    if(is_cubic){
        for(size_t i = 0; i < N_; i++){
            for(size_t j = 0; j < 3; j++){
                A_cubic_[i*3 + j] = G[i][j];
            }
        }
    }else{
        for(size_t i = 0; i < Nf_; i++){
            deg_[i] = G[i].size();
            for(size_t j = 0; j < G[i].size(); j++){
                A_dual_[i*6 + j] = G[i][j];
            }
        }
    }
    if(allocate_hessian_etc){
        hessian_ = SyclVector<T>(N_*90);
        eigenvalues_ = SyclVector<T>(N_*3);
        eigenvectors_ = SyclVector<T>(N_*3*N_*3);
    }
}

template <typename T, typename K>
FullereneIsomer<T,K>::FullereneIsomer(const Graph& G, bool is_cubic, bool allocate_hessian_etc) : 
    FullereneIsomer(G.neighbours, is_cubic, allocate_hessian_etc) {}

template <typename T, typename K>
FullereneIsomer<T,K>::FullereneIsomer(const Polyhedron& P, bool is_cubic, bool allocate_hessian_etc) : 
    FullereneIsomer(P.neighbours, is_cubic, allocate_hessian_etc) {
    if(is_cubic){
        for(size_t i = 0; i < N_; i++){
            for(size_t j = 0; j < 3; j++){
                X_cubic_[i*3 + j] = P.points[i][j];
            }
        }
    } else {
        for(size_t i = 0; i < Nf_; i++){
            for(size_t j = 0; j < 3; j++){
                X_dual_[i*3 + j] = P.points[i][j];
            }
        }
    }
}

template <typename T, typename K>
FullereneIsomer<T,K>::FullereneIsomer(const FullereneIsomer<T,K>& other) = default;

template <typename T, typename K>
FullereneIsomer<T,K>::FullereneIsomer(FullereneIsomer<T,K>&& other) = default;

template <typename T, typename K>
FullereneIsomer<T,K>& FullereneIsomer<T,K>::operator=(const FullereneIsomer<T,K>& other) = default;




DeviceWrapper::DeviceWrapper(size_t id, DeviceType device_type) : id_(id), type_(device_type) {
    switch(type_){
        case DeviceType::CPU:
            device_ = static_cast<const void*>(reinterpret_cast<sycl::device*>(__FULLERENE_ALLCPUDEVICES__.data()) + id);
            break;
        case DeviceType::GPU:
            device_ = static_cast<const void*>(reinterpret_cast<sycl::device*>(__FULLERENE_ALLGPUDEVICES__.data()) + id);
            break;
        case DeviceType::ACCELERATOR:
            device_ = static_cast<const void*>(reinterpret_cast<sycl::device*>(__FULLERENE_ALLACCELERATORS__.data()) + id);
            break;
    }
}
    
bool DeviceWrapper::is_cpu() const {return type_ == DeviceType::CPU;}
bool DeviceWrapper::is_gpu() const {return type_ == DeviceType::GPU;}

size_t DeviceWrapper::get_id() const  {return id_;}  
std::string DeviceWrapper::get_name() const  {return DEVICE_CAST(device_).get_info<sycl::info::device::name>();}
    
size_t DeviceWrapper::get_property(DeviceProperty property) const {
    static std::array<size_t, (int)DeviceProperty::NUMBER_OF_PROPERTIES> properties{ //These must be in the same order as DeviceProperty enum
        DEVICE_CAST(device_).get_info<sycl::info::device::max_work_group_size>(),
        DEVICE_CAST(device_).get_info<sycl::info::device::max_clock_frequency>(),
        DEVICE_CAST(device_).get_info<sycl::info::device::max_compute_units>(),
        DEVICE_CAST(device_).get_info<sycl::info::device::max_mem_alloc_size>(),
        DEVICE_CAST(device_).get_info<sycl::info::device::global_mem_size>(),
        DEVICE_CAST(device_).get_info<sycl::info::device::local_mem_size>(),
        DEVICE_CAST(device_).get_info<sycl::info::device::max_constant_buffer_size>(),
        DEVICE_CAST(device_).get_info<sycl::info::device::max_constant_args>()
        };
    return properties[(int)property];}


QueueWrapper::QueueWrapper(int device_id, const DeviceType device_type, bool in_order) : device_(device_id, device_type){
    const sycl::device& device = device_type == DeviceType::CPU ? reinterpret_cast<sycl::device*>(DeviceWrapper::__FULLERENE_ALLCPUDEVICES__.data())[device_id] : reinterpret_cast<sycl::device*>(DeviceWrapper::__FULLERENE_ALLGPUDEVICES__.data())[device_id];
    auto tempqueue = in_order ? sycl::queue(device, sycl::property::queue::in_order{}) : sycl::queue(device);
    queue_.resize(sizeof(sycl::queue));
    memcpy((char*)queue_.data(), (char*)&tempqueue, sizeof(sycl::queue));
    id_ = device_id;
}
void QueueWrapper::wait() {reinterpret_cast<sycl::queue*>(queue_.data())->wait();}

const DeviceWrapper QueueWrapper::get_device() const { return device_;}

template struct SyclVector<float>;
template struct SyclVector<double>;
template struct SyclVector<uint16_t>;
template struct FullereneIsomer<float,uint16_t>;
template struct FullereneIsomer<double,uint16_t>;

