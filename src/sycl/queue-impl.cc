#pragma once

#include <CL/sycl.hpp>
#include <fullerenes/sycl-wrappers.hh>
#ifndef DEVICE_CAST
    #define DEVICE_CAST(x,ix) (reinterpret_cast<const sycl::device*>(x)[ix])
#endif

struct SyclQueueImpl : public sycl::queue{
    using sycl::queue::queue;

    inline static const auto device_arrays = std::array{ 
                sycl::device::get_devices(sycl::info::device_type::cpu), 
                sycl::device::get_devices(sycl::info::device_type::gpu), 
                sycl::device::get_devices(sycl::info::device_type::accelerator) };

    SyclQueueImpl(Device dev, bool in_order) : sycl::queue( device_arrays.at((int)dev.second).at(dev.first), in_order ? sycl::property::queue::in_order{} : sycl::property_list{}) {}

    operator SyclQueue() {return SyclQueue(std::make_unique<SyclQueueImpl>(std::move(*this)));}
};

struct DeviceImpl : public sycl::device{
    using sycl::device::device;
    DeviceImpl(Device dev) : sycl::device(SyclQueueImpl::device_arrays.at((int)dev.second).at(dev.first)) {}
};