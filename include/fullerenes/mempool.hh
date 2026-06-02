#pragma once
#include <cstdint>
#include <stdexcept>
#include <fullerenes/sycl-headers/sycl-device-queue.hh>
#include <fullerenes/sycl-headers/sycl-span.hh>

struct BumpAllocator {
    template <typename T>
    BumpAllocator(T* data, size_t byte_size): data(data), byte_size(byte_size){}

    template <typename T>
    BumpAllocator(std::span<T> span): data(span.data()), byte_size(span.size()*sizeof(T)){}

    template<typename T>
    constexpr inline static auto alignment(const Device& device){
        auto device_align_bytes = device.get_property(DeviceProperty::MEM_BASE_ADDR_ALIGN) / 8;
        return std::max(device_align_bytes, static_cast<std::uintptr_t>(alignof(T)));
    }

    template<typename T>
    constexpr inline static size_t allocation_size(const Device& device, size_t size){
        std::uintptr_t total_size = size * sizeof(T);
        return (total_size + alignment<T>(device) - 1) & ~(alignment<T>(device) - 1);
    }

    template<typename T>
    constexpr inline static size_t allocation_size(SyclQueue& ctx, size_t size)   {return allocation_size<T>(ctx.device(), size);}

    template<typename T>
    constexpr inline std::span<T> allocate(const Device& device, size_t size){
        size_t alloc_size = allocation_size<T>(device,size);
        if (alloc_size > byte_size){
            throw std::runtime_error("Out of memory");
        }
        
        std::uintptr_t addr = reinterpret_cast<std::uintptr_t>(data);
        std::uintptr_t aligned = (addr + alignment<T>(device) - 1) & ~(alignment<T>(device) - 1);  // Changed from alignof(T)
        T* ptr = reinterpret_cast<T*>(aligned);

        data = reinterpret_cast<void*>(ptr + size);
        byte_size -= (reinterpret_cast<char*>(data) - reinterpret_cast<char*>(ptr));

        return std::span(ptr, size);
    }

    template<typename T>
    constexpr inline std::span<T> allocate(SyclQueue& ctx, size_t size) {return allocate<T>(ctx.device(), size);}

    // Bytes consumed since construction (for debug-asserting that
    // *_buffer_size() queries match what the kernel actually allocates).
    constexpr inline size_t bytes_used() const { return total_byte_size - byte_size; }
    constexpr inline size_t bytes_remaining() const { return byte_size; }

    private:

        void*  data;
        size_t byte_size;
        size_t total_byte_size = byte_size;
};