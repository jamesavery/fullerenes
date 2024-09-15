#pragma once
#include <string>
#include <cstdint>
#include <array>
#include <span>
#include <memory>
#include <execution>
#include <numeric>
#include <fullerenes/graph.hh>
#include <fullerenes/polyhedron.hh>
#include <fullerenes/sycl-headers/sycl-device-queue.hh>
#include <fullerenes/sycl-headers/sycl-span.hh>
#include <fullerenes/sycl-headers/sycl-vector.hh>
#include <fullerenes/sycl-headers/sycl-fullerene-structs.hh>

using uint64_t = unsigned long;
using uint32_t = unsigned int;
using uint16_t = unsigned short;
using uint8_t  = unsigned char;

//Use this definition in class definitions to add the following data members

template <typename Func, typename... Args>
auto execute_with_policy(SyclQueue& Q, const Policy policy, Func&& func, Args&&... args){
    if(policy == Policy::SYNC){
        func(std::forward<Args>(args)...);
        Q.wait();
    } else {
        func(std::forward<Args>(args)...);
    }
}