#pragma once

#include <CL/sycl.hpp>
#include <fullerenes/sycl-headers/sycl-device-queue.hh>

#ifndef DEVICE_CAST
    #define DEVICE_CAST(x,ix) (reinterpret_cast<const sycl::device*>(x)[ix])
#endif

struct SyclEventImpl : public sycl::event{
    using sycl::event::event;

    SyclEventImpl(sycl::event&& event) : sycl::event(event) {}
};

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



struct SyclQueueImpl : public sycl::queue{
    using sycl::queue::queue;

    inline static const auto device_arrays = std::array{ 
                sycl::device::get_devices(sycl::info::device_type::cpu), 
                sycl::device::get_devices(sycl::info::device_type::gpu), 
                sycl::device::get_devices(sycl::info::device_type::accelerator),
                sycl::device::get_devices(sycl::info::device_type::host)};

    static_assert(device_arrays.size() == (int)DeviceType::NUM_DEV_TYPES && "DeviceType enum does not match device_arrays size");
    
    SyclQueueImpl(Device dev, bool in_order) : sycl::queue( device_arrays.at((int)dev.type).at(dev.idx), in_order ? sycl::property::queue::in_order{} : sycl::property_list{}), device_(dev) {}
    SyclQueueImpl() : sycl::queue( device_arrays.at((int)DeviceType::CPU).at(0)), device_(Device{0, DeviceType::CPU}) {}
    const Device device_;
};

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
        h.depends_on(static_cast<sycl::event>(*event));
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
        default: std::cerr << "Unknown property" << std::endl; return 0;
    }
}
