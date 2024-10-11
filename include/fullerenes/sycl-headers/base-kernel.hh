#pragma once
//#include <fullerenes/sycl-isomer-batch.hh>
#include <fullerenes/sycl-headers/sycl-device-queue.hh>
#include <fullerenes/sycl-headers/sycl-fullerene-structs.hh>
#include <fullerenes/sycl-headers/sycl-vector.hh>
#include <fullerenes/sycl-headers/sycl-span.hh>
#include <fullerenes/sycl-headers/misc-enums.hh>
#include <thread>
#include <future>

template <template <typename...> class Template, typename T>
struct is_specialization_of : std::false_type {};

template <template <typename...> class Template, typename... Args>
struct is_specialization_of<Template, Template<Args...>> : std::true_type {};

template <typename T>
auto partition_vector(size_t ix, size_t batch_size, const T& value) {
    if constexpr (is_specialization_of<SyclVector, T>::value || is_specialization_of<Span, T>::value) {
        auto size = value.size();
        if (size % batch_size != 0) {
            throw std::runtime_error("The size of the SyclVector is not divisible by the batch size");
        }
        auto size_per_isomer = size / batch_size;
        auto modified = Span(value.begin() + ix*size_per_isomer, value.begin() + (ix+1)*size_per_isomer);
        return modified;
    } else {
        return value;  // Leave other types unchanged
    }
}

// Function to modify the parameter pack
template <typename... Args>
auto prepare_execution_arguments(size_t ix, size_t batch_size, const Args&... args) {
    return std::make_tuple(partition_vector(ix, batch_size, args)...);  // Apply the modification to each argument
}

template <typename T>
using ArrayOfSyclVectors = std::array<std::vector<std::vector<SyclVector<T>>>, (size_t)DeviceType::NUM_DEV_TYPES>;


template <typename T>
struct DefaultInitAtomic : std::atomic<T> {
    using std::atomic<T>::atomic;
    DefaultInitAtomic() : std::atomic<T>(T{}) {}
    DefaultInitAtomic(const DefaultInitAtomic<T>& other) : std::atomic<T>(other.load()) {}
    DefaultInitAtomic(T desired) : std::atomic<T>(desired) {}
};



struct DispatchCounters : public std::array<std::vector<DefaultInitAtomic<int>>, (size_t)DeviceType::NUM_DEV_TYPES> {
    DispatchCounters() = default;
    DispatchCounters(const DispatchCounters& other) = default;
    DispatchCounters(DispatchCounters&& other) = default;
    DispatchCounters& operator=(const DispatchCounters& other) = default;

    auto& operator[](Device device) {
        return this->at((size_t)device.type).at(device.idx);
    }    
    auto& operator[](SyclQueue& Q) {
        return this->at((size_t)Q.device_.type).at(Q.device_.idx);
    }
};

template <typename T>
struct FunctorArrays : public ArrayOfSyclVectors<T> {
    auto& operator[](const Device device) {
        return this->at((size_t)device.type).at(device.idx);
    }
    auto& operator[](const SyclQueue& Q) {
        return (*this)[Q.device_];
    }

    auto& operator[](const DeviceType device_type) {
        return this->at((size_t)device_type);
    }

    Span<T> operator[](const std::pair<Device, size_t> ix_tuple) {
        return this->at((size_t)ix_tuple.first.type).at(ix_tuple.first.idx).at(ix_tuple.second);
    }
    Span<T> operator[](const std::pair<SyclQueue&, size_t> ix_tuple) {
        auto& [Q, ix] = ix_tuple;
        return (*this)[std::pair<Device,size_t>{Q.device_, ix}];
    }
};

struct MutexEvent : public std::tuple<std::mutex, SyclEvent> {
    using std::tuple<std::mutex, SyclEvent>::tuple;

    void unlock_and_wait() {
        std::get<std::mutex>(*this).unlock();
        std::get<SyclEvent>(*this).wait();
    }

    void lock_and_wait() {
        std::get<std::mutex>(*this).lock();
        std::get<SyclEvent>(*this).wait();
    }

    void unlock() {
        std::get<std::mutex>(*this).unlock();
    }

    void wait() {
        std::get<SyclEvent>(*this).wait();
    }

    template <typename EventT>
    void set(EventT&& event) {
        std::get<SyclEvent>(*this) = event;
    }

    MutexEvent() = default;
    MutexEvent(const MutexEvent& other) = default;
    MutexEvent(MutexEvent&& other) = default;
    MutexEvent& operator=(const MutexEvent& other) = default;
    MutexEvent& operator=(MutexEvent&& other) = default;

    MutexEvent(SyclEvent&& event)  {
        std::get<SyclEvent>(*this) = std::move(event);
    }

    MutexEvent& operator=(SyclEvent&& event) {
        std::get<SyclEvent>(*this) = std::move(event);
        return *this;
    }
};

struct MutexVector : public std::vector<MutexEvent> {
    void wait() {
        for (auto& mut : *this) {
            mut.wait();
        }
    }

    void unlock_and_wait() {
        for (auto& mut : *this) {
            mut.unlock_and_wait();
        }
    }

    void lock_and_wait() {
        for (auto& mut : *this) {
            mut.lock_and_wait();
        }
    }

    void unlock() {
        for (auto& mut : *this) {
            mut.unlock();
        }
    }


    void resize(size_t N) {
        this->unlock_and_wait();
        this->clear();
        auto new_vector = std::vector<MutexEvent>(N);
        this->swap(new_vector);
    }
};

template <typename T>
struct AOVOV : public std::array<std::vector<std::vector<T>>, (size_t)DeviceType::NUM_DEV_TYPES> {
    auto& operator[](std::pair<Device, size_t> ix_tuple) {
        return this->at((size_t)ix_tuple.first.type).at(ix_tuple.first.idx).at(ix_tuple.second);
    }    
    auto& operator[](std::pair<SyclQueue&, size_t> ix_tuple) {
        auto& [Q, ix] = ix_tuple;
        return this->at((size_t)Q.device_.type).at(Q.device_.idx).at(ix);
    }

    auto& operator[](SyclQueue& Q) {
        return this->at((size_t)Q.device_.type).at(Q.device_.idx);
    }

    auto& operator[](Device device) {
        return this->at((size_t)device.type).at(device.idx);
    }

    auto& operator[](DeviceType device_type) {
        return this->at((size_t)device_type);
    }

    auto& operator[](int device_type) {
        return (*this)[static_cast<DeviceType>(device_type)];
    }
};

struct FunctorMutexes : public std::array<std::vector<MutexVector>, (size_t)DeviceType::NUM_DEV_TYPES> {
    auto& operator[](std::pair<Device, size_t> ix_tuple) {
        return this->at((size_t)ix_tuple.first.type).at(ix_tuple.first.idx).at(ix_tuple.second);
    }    
    auto& operator[](std::pair<SyclQueue&, size_t> ix_tuple) {
        auto& [Q, ix] = ix_tuple;
        return this->at((size_t)Q.device_.type).at(Q.device_.idx).at(ix);
    }

    auto& operator[](SyclQueue& Q) {
        return this->at((size_t)Q.device_.type).at(Q.device_.idx);
    }

    auto& operator[](Device device) {
        return this->at((size_t)device.type).at(device.idx);
    }

    auto& operator[](DeviceType device_type) {
        return this->at((size_t)device_type);
    }

    auto& operator[](int device_type) {
        return (*this)[static_cast<DeviceType>(device_type)];
    }

    void wait() {
        for (auto& dev_type : *this) {
            for (auto& device_events : dev_type) {
                for (auto& mut : device_events) {
                    mut.wait();
                }
            }
        }
    }
};

template <typename KernelImpl>
struct KernelFunctor{
    DispatchCounters dispatch_counters_;
    FunctorMutexes mutexes_;
    AOVOV<std::future<SyclEvent>> futures_;
    bool initialized = false;

    KernelFunctor() {
        for (int i = 0; i < dispatch_counters_.size(); i++) {
            dispatch_counters_.at(i).resize(Device::get_devices(static_cast<DeviceType>(i)).size());
            mutexes_.at(i).resize(Device::get_devices(static_cast<DeviceType>(i)).size());
            futures_.at(i).resize(Device::get_devices(static_cast<DeviceType>(i)).size());
        }
    }


    //Default implementation of get_max_concurrent_launches, can be overridden by each functor.
    template <typename... Args>
    inline constexpr size_t get_max_concurrent_launches(SyclQueue& Q, size_t N, Args&&... args) const {
        size_t max_cus = Q.device_.get_property(DeviceProperty::MAX_COMPUTE_UNITS);
        size_t max_wg_size = Q.device_.get_property(DeviceProperty::MAX_WORK_GROUP_SIZE);
        size_t num_wgs_per_N = (N + max_wg_size - 1) / max_wg_size;
        if (num_wgs_per_N > max_cus) return 1;
        return max_cus / num_wgs_per_N;
    }

    template <typename... Args>
    SyclEvent compute(Args&&... args) const {
        throw std::logic_error("KernelFunctor::compute() not implemented");
    }

    template <typename... Args>
    inline constexpr auto to_tuple(size_t N, Args&&... args) const {
        return static_cast<const KernelImpl*>(this)->to_tuple(N, std::forward<Args>(args)...);
    }

    template <typename... Args>
    inline constexpr auto to_tuple_batch(size_t N, Args&&... args) const {
        return static_cast<const KernelImpl*>(this)->to_tuple_batch(N, std::forward<Args>(args)...);
    }

    template <typename... Args>
    void initialize(Args&&... args) {
        if (initialized) return;
        auto pair_tuple = to_tuple(0, std::forward<Args>(args)...);
        std::apply([&](auto&... pairs) {
            auto init_and_check_capacity = [&](auto& array_of_vectors) {
                for (int i = 0; i < array_of_vectors.size(); i++) {
                    auto n_device_of_type = Device::get_devices(static_cast<DeviceType>(i)).size();
                    (array_of_vectors[static_cast<DeviceType>(i)]).resize(n_device_of_type);
                }
            };
            (init_and_check_capacity(pairs.first), ...);
        }, pair_tuple);
        initialized = true;
    }
    
    template <bool is_batch, typename... Args>
    inline constexpr auto allocate_and_return_tuple(SyclQueue& Q, size_t launch_idx, size_t N, size_t batch_size, Args&&... args) {
        initialize(std::forward<Args>(args)...);
        using TupleType = typename std::conditional_t<is_batch, decltype(to_tuple_batch(N, std::forward<Args>(args)...)), decltype(to_tuple(N, std::forward<Args>(args)...))>;
        TupleType pair_tuple = [this, N](auto&&... args) -> TupleType {
            if constexpr (is_batch) {
                return to_tuple_batch(N, args...);
            }else{
                return to_tuple(N, args...);
            }
        }(std::forward<Args>(args)...);

        std::apply([&](auto&... pairs) {
            auto allocate_and_check_capacity = [&](auto& array_of_vectors, size_t capacity) {
                array_of_vectors[Q].resize(get_max_concurrent_launches(Q, N, std::forward<Args>(args)...));
                array_of_vectors[Q].at(launch_idx).resize(capacity*batch_size);
            };
            (allocate_and_check_capacity(pairs.first, pairs.second), ...);
        }, pair_tuple);
        return std::apply([&Q, launch_idx](auto&... pairs) {
            return std::make_tuple(pairs.first[{Q, launch_idx}]...);
        }, pair_tuple);
    }

    template <typename T, typename K, typename... Args>
    auto operator() (SyclQueue& Q, FullereneBatchView<T,K> batch, LaunchPolicy policy, Args&&... args) {
        if (policy == LaunchPolicy::SYNC) Q.wait();
        auto max_concurrent_launches = static_cast<KernelImpl*>(this)->get_max_concurrent_launches(Q, batch.N_, std::forward<Args>(args)...);
        mutexes_[Q].resize(std::min((size_t)batch.size(), max_concurrent_launches));
        futures_[Q].resize(std::min((size_t)batch.size(), max_concurrent_launches));

        if (batch.size() == 0) {return;}
        if (batch.N_ > Q.device_.get_property(DeviceProperty::MAX_WORK_GROUP_SIZE)){
            SyclQueue out_of_order_queue(Q.device_, false);
            auto full_ix = 0;
            std::for_each(batch.begin(), batch.end(), [&](auto isomer) {
                auto circular_ix = (this->dispatch_counters_[out_of_order_queue]++) % std::min((size_t)batch.size(), max_concurrent_launches);
                //Mutex ensures that no two try to resize the same vector at the same time, nor use the same memory for different isomers
                //mutexes_[{out_of_order_queue, circular_ix}].lock_and_wait();
                //if (threads_[{out_of_order_queue, circular_ix}].joinable()) threads_[{out_of_order_queue, circular_ix}].join();
                /* threads_[{out_of_order_queue, circular_ix}] = std::thread([=](SyclQueue& lambda_queue, KernelFunctor& kernel){
                    auto data_tuple = kernel.allocate_and_return_tuple(lambda_queue, circular_ix, isomer.N_);
                    isomer_function(out_of_order_queue, isomer, data_tuple);
                }, std::ref(out_of_order_queue), std::ref(this)); */
                
                auto data_tuple = allocate_and_return_tuple<false>(out_of_order_queue, circular_ix, isomer.N_, 1, std::forward<Args>(args)...);
                auto prepared_args = prepare_execution_arguments(full_ix, batch.size(), std::forward<Args>(args)...);
                std::apply([&](auto&&... args) {
                    static_cast<KernelImpl*>(this)->compute(out_of_order_queue, isomer, std::forward<decltype(args)>(args)...);
                }, std::tuple_cat(prepared_args, data_tuple));
                //isomer_function(out_of_order_queue, isomer, std::tuple_cat(data_tuple, prepared_args));
                //auto data_tuple = allocate_and_return_tuple(out_of_order_queue, circular_ix, batch.N_);
                //mutexes_[{out_of_order_queue, circular_ix}] = isomer_function(out_of_order_queue, isomer, data_tuple);
                //mutexes_[{out_of_order_queue, circular_ix}].unlock();
                full_ix++;
            });
            //Enqueue all the events in the input queue, this way asynchronicity is preserved
            //When the caller waits for the input queue, it will wait for all the out_of_order_queue events
            /* std::for_each(futures_[Q].begin(), futures_[Q].end(), [&](auto& fut) {
                fut.wait();
            }); */
        }else{
            auto batch_data = allocate_and_return_tuple<true>(Q, 0, batch.N_, batch.size(), std::forward<Args>(args)...);
            //batch_function(Q, batch, batch_data);
            std::apply([&](auto&&... data) {
                static_cast<KernelImpl*>(this)->compute(Q, batch, std::forward<Args>(args)..., std::forward<decltype(data)>(data)...);
            }, batch_data);
        }
        if (policy == LaunchPolicy::SYNC) Q.wait();
    }

    template <typename T, typename K, typename... Args>
    auto operator() (SyclQueue& Q, FullereneBatch<T,K>& batch, LaunchPolicy policy, Args&&... args) {
        this->operator()(Q, (FullereneBatchView<T,K>)batch, policy, std::forward<Args>(args)...);
    }

    template <typename T, typename K, typename... Args>
    auto operator() (SyclQueue& Q, Fullerene<T,K> isomer, LaunchPolicy policy, Args&&... args) {
        if (policy == LaunchPolicy::SYNC) Q.wait();
        auto ret_val = std::apply([&](auto&&... data) {
            return static_cast<KernelImpl*>(this)->compute(Q, isomer, std::forward<Args>(args)..., data...);
        }, allocate_and_return_tuple<false>(Q, 0, isomer.N_, 1, std::forward<Args>(args)...));
        if (policy == LaunchPolicy::SYNC) Q.wait();
        return ret_val;
    }
};