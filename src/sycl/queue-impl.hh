#pragma once

#include <sycl/sycl.hpp>
#include <fullerenes/sycl-headers/sycl-device-queue.hh>
#include <tuple>
#include <type_traits>

#ifndef DEVICE_CAST
    #define DEVICE_CAST(x,ix) (reinterpret_cast<const sycl::device*>(x)[ix])
#endif

struct SyclEventImpl {
    sycl::event event;

    SyclEventImpl() = default;
    SyclEventImpl(sycl::event event_) : event(std::move(event_)) {}

    void wait() const { const_cast<sycl::event&>(event).wait(); }

    operator sycl::event&() { return event; }
    operator const sycl::event&() const { return event; }
};

#ifdef DEFINE_SYCL_QUEUE_METHODS
SyclEvent::SyclEvent() : impl_(std::make_unique<SyclEventImpl>(sycl::event())) {}
SyclEvent::SyclEvent(SyclEventImpl&& impl) : impl_(std::make_unique<SyclEventImpl>(std::move(impl))) {}
SyclEvent& SyclEvent::operator=(SyclEvent&& other) {
    impl_ = std::move(other.impl_);
    return *this;
}
SyclEvent& SyclEvent::operator=(SyclEventImpl&& impl) {
    impl_ = std::make_unique<SyclEventImpl>(std::move(impl));
    return *this;
}

SyclEvent::SyclEvent(SyclEvent&& other) = default;


SyclEvent::~SyclEvent() = default;

void SyclEvent::wait() const {impl_->wait();}

SyclEventImpl* SyclEvent::operator ->() const {return impl_.get();}
SyclEventImpl& SyclEvent::operator *() const {return *impl_;}

#endif // DEFINE_SYCL_QUEUE_METHODS

struct SyclQueueImpl : public sycl::queue{
    using sycl::queue::queue;

    // Heap-allocated and intentionally never freed: the cached sycl::device
    // objects keep the AdaptiveCpp runtime alive, and destroying them from a
    // static destructor at process exit tears the runtime down out of order
    // (SIGSEGV in allocation_tracker::unregister_allocation during
    // __cxa_finalize). Leaking this one-time cache lets the OS reclaim it at
    // process death without the crashing teardown.
    inline static const auto& device_arrays = *new std::array{ 
                sycl::device::get_devices(sycl::info::device_type::cpu), 
                sycl::device::get_devices(sycl::info::device_type::gpu), 
                sycl::device::get_devices(sycl::info::device_type::accelerator),
                sycl::device::get_devices(sycl::info::device_type::host)};

    // size() on the runtime-initialized reference is not constexpr; check the
    // (CTAD-deduced) array size from its type instead.
    static_assert(std::tuple_size_v<std::remove_reference_t<decltype(device_arrays)>> == (size_t)DeviceType::NUM_DEV_TYPES && "DeviceType enum does not match device_arrays size");
    
    SyclQueueImpl(Device dev, bool in_order) : sycl::queue( device_arrays.at((int)dev.type).at(dev.idx), in_order ? sycl::property::queue::in_order{} : sycl::property_list{}), device_(dev) {}
    SyclQueueImpl() : sycl::queue( device_arrays.at((int)DeviceType::CPU).at(0)), device_(Device{0, DeviceType::CPU}) {}
    const Device device_;
};

// @anchor sycl-launch-per-isomer
// Launch one work-group per isomer over a batch: the work-group size is
// work_per_isomer and there are `count` isomers, so the global range is
// work_per_isomer*count and each group handles one isomer. Owns the empty-batch
// guard, the nd_range, the submit, and the SyclEvent wrap -- the boilerplate every
// batch kernel shares. The command-group function cgf(handler&, nd_range<1>) creates
// the kernel's local accessors and issues h.parallel_for(ndr, ...).
// @post on empty:    (count == 0 || work_per_isomer == 0) -> an already-complete event,
//                        no kernel submitted. A zero global range is a no-op on the
//                        host/OpenMP backend but cuLaunchKernel rejects a zero grid
//                        with CUDA_ERROR_INVALID_VALUE, so empty batches must not launch.
// @post on nonempty: submits cgf(h, nd_range(work_per_isomer*count, work_per_isomer))
//                        and returns the event for that submission.
template <class CGF>
inline SyclEvent launch_per_isomer(SyclQueue& Q, size_t work_per_isomer, size_t count, CGF&& cgf) {
    if (count == 0 || work_per_isomer == 0) return SyclEvent();
    sycl::nd_range<1> ndr(sycl::range<1>(work_per_isomer * count), sycl::range<1>(work_per_isomer));
    SyclEventImpl ev = Q->submit([&](sycl::handler& h){ cgf(h, ndr); });
    return SyclEvent(std::move(ev));
}

#ifdef DEFINE_SYCL_QUEUE_METHODS
SyclQueue::SyclQueue() : device_({0, DeviceType::CPU}), in_order_(true) {
    impl_ = std::make_unique<SyclQueueImpl>(device_, in_order_);
}

SyclQueue::SyclQueue(Device dev, bool in_order) : device_(dev), in_order_(in_order) {
    impl_ = std::make_unique<SyclQueueImpl>(dev, in_order);
}

SyclQueue::~SyclQueue() = default;
SyclQueue::SyclQueue(SyclQueue&& other) = default;
SyclQueue& SyclQueue::operator=(SyclQueue&& other) = default;

void SyclQueue::wait() const {impl_->wait();}
void SyclQueue::wait_and_throw() const {impl_->wait_and_throw();}


SyclQueueImpl* SyclQueue::operator->() const {
    return impl_.get();
}

SyclQueueImpl& SyclQueue::operator*() const {
    return *impl_;
}

void SyclQueue::enqueue(SyclEvent& event) {
    if(event.impl_.get() == nullptr) return;
    impl_->submit([&](sycl::handler& h) { 
        h.depends_on(event.impl_->event);
    });
}

SyclEvent SyclQueue::get_event() {
    SyclEventImpl event = impl_->submit([](sycl::handler& h){h.single_task([](){});});
    return event;
}

std::vector<Device> Device::get_devices(DeviceType type){
    std::vector<Device> devices(SyclQueueImpl::device_arrays.at(static_cast<int>(type)).size());
    std::generate(devices.begin(), devices.end(), 
        [i = 0, type]() mutable { return Device(i,type); });
    return devices;
}


std::string Device::get_name() const  {
    return SyclQueueImpl::device_arrays.at(static_cast<int>(type)).at(idx).get_info<sycl::info::device::name>();
}

std::string Device::get_vendor() const {
    return SyclQueueImpl::device_arrays.at(static_cast<int>(type)).at(idx).get_info<sycl::info::device::vendor>();
}

size_t Device::get_property(DeviceProperty property) const {
    auto d = SyclQueueImpl::device_arrays.at(static_cast<int>(type)).at(idx);
    
    switch(int(property)){
        case 0: return d.get_info<sycl::info::device::max_work_group_size>();
        case 1: return d.get_info<sycl::info::device::max_clock_frequency>();
        case 2: return d.get_info<sycl::info::device::max_compute_units>();
        case 3: return d.get_info<sycl::info::device::max_mem_alloc_size>();
        case 4: return d.get_info<sycl::info::device::global_mem_size>();
        case 5: return d.get_info<sycl::info::device::local_mem_size>();
        case 7: return d.get_info<sycl::info::device::max_constant_args>();
        case 8: return d.get_info<sycl::info::device::max_num_sub_groups>();
        case 9: return d.get_info<sycl::info::device::sub_group_sizes>()[0];
        case 10: return d.get_info<sycl::info::device::mem_base_addr_align>();
        default: std::cerr << "Unknown property" << std::endl; return 0;
    }
}
#endif // DEFINE_SYCL_QUEUE_METHODS
