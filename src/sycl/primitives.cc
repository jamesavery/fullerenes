
#ifdef __INTEL_LLVM_COMPILER_FIX_THE_ONEAPI_DPL_INCLUDES_ // When the oneDPL is fixed, we can change this to "__INTEL_LLVM_COMPILER"
    #include <oneapi/dpl/algorithm>
    #include <oneapi/dpl/execution>
    #include <oneapi/dpl/iterator>
    #include <oneapi/dpl/numeric>
#endif
#include "primitives.hh"
#include <fullerenes/sycl-wrappers.hh>
    


#ifdef __INTEL_LLVM_COMPILER_FIX_THE_ONEAPI_DPL_INCLUDES_
#define DEFINE_FUNCS(templates, returntype, funcname, post_args) \
    templates returntype funcname(SyclQueue& queue, T* begin, T* end, T* result, post_args){\
        auto& Q = queue.impl_->get_queue();\
        auto policy = oneapi::dpl::execution::make_device_policy(Q.get_device());\
        return oneapi::dpl::funcname(policy, begin, end, result, post_args);\
    }\
    templates returntype funcname(SyclQueue& queue, const SyclVector<T>& vec, SyclVector<T>& result, post_args){\
        assert(vec.capacity() == result.capacity());\
        return funcname(queue, vec.data(), vec.data() + vec.capacity(), result.data(), post_args);\
    }\
    templates returntype funcname(SyclQueue& queue, const Span<T> vec, const Span<T> result, post_args){\
        assert(vec.size() == result.size());\
        return funcname(queue, vec.data(), vec.data() + vec.size(), result.data(), post_args);\
    }\
    templates returntype funcname(SyclQueue& queue, SyclVector<T>& vec, post_args){\
        return funcname(queue, vec.data(), vec.data() + vec.size(), vec.data(), post_args);\
    }\
    templates returntype funcname(SyclQueue& queue, Span<T> vec, post_args){\
        return funcname(queue, vec.data(), vec.data() + vec.size(), vec.data(), post_args);\
    }
#else
#define DEFINE_FUNCS(templates, returntype, funcname, post_args, call_args) \
    templates returntype funcname(SyclQueue& queue, T* begin, T* end, T* result, post_args){\
        if constexpr (std::is_void_v<returntype>){\
            std::funcname(std::execution::par_unseq, begin, end, result, call_args);\
        }else{\
            return std::funcname(std::execution::par_unseq, begin, end, result, call_args);}\
    }\
    templates returntype funcname(SyclQueue& queue, const SyclVector<T>& vec, SyclVector<T>& result, post_args){\
        assert(vec.capacity() == result.capacity());\
        if constexpr (std::is_void_v<returntype>){\
            funcname(queue, vec.data(), vec.data() + vec.capacity(), result.data(), call_args);\
        }else{\
            return funcname(queue, vec.data(), vec.data() + vec.capacity(), result.data(), call_args);}\
    }\
    templates returntype funcname(SyclQueue& queue, const Span<T> vec, const Span<T> result, post_args){\
        assert(vec.size() == result.size());\
        if constexpr (std::is_void_v<returntype>){\
            funcname(queue, vec.data(), vec.data() + vec.size(), result.data(), call_args);\
        }else{\
            return funcname(queue, vec.data(), vec.data() + vec.size(), result.data(), call_args);}\
    }\
    templates returntype funcname(SyclQueue& queue, SyclVector<T>& vec, post_args){\
        if constexpr (std::is_void_v<returntype>){\
            funcname(queue, vec.data(), vec.data() + vec.size(), vec.data(), call_args);\
        }else{\
            return funcname(queue, vec.data(), vec.data() + vec.size(), vec.data(), call_args);}\
    }\
    templates returntype funcname(SyclQueue& queue, Span<T> vec, post_args){\
        if constexpr (std::is_void_v<returntype>){\
            funcname(queue, vec.data(), vec.data() + vec.size(), vec.data(), call_args);\
        }else{\
            return funcname(queue, vec.data(), vec.data() + vec.size(), vec.data(), call_args);}\
    }
#endif





namespace primitives{
    DECLARE_OR_DEFINE_ALL_FUNCS(DEFINE_FUNCS)

    template <typename B>
    constexpr void my_assert() { 
        static_assert(false, "oh no");
    }

    template<typename Element1, typename Element2>
    void iterateWithTwoElements(Element1, Element2, std::tuple<>) {}
    template<typename Element1, typename Element2, typename... Pack3>
    auto iterateWithTwoElements(Element1, Element2, std::tuple<Pack3...> pack3) {
        return std::make_tuple( 
        [](auto... args) 
        { 
            return exclusive_scan<Element1, Element2, Pack3>(args...); 
        }...);
    }
    template<typename Element1, typename... Pack2, typename... Pack3>
    void iterateWithElement(Element1, std::tuple<Pack2...> pack2, std::tuple<Pack3...> pack3) {
        ([&](auto& elem2){
            iterateWithTwoElements(Element1{}, elem2, pack3);
        }(std::get<Pack2>(pack2)), ...);
    }
    template<typename... Pack2, typename... Pack3>
    void iteratePacks(std::tuple<>, std::tuple<Pack2...>, std::tuple<Pack3...>) {}

    template<typename Head1, typename... Tail1, typename... Pack2, typename... Pack3>
    void iteratePacks(std::tuple<Head1, Tail1...> pack1, std::tuple<Pack2...> pack2, std::tuple<Pack3...> pack3) {
        iterateWithElement(std::get<0>(pack1), pack2, pack3);
        iteratePacks(std::tuple<Tail1...>{}, pack2, pack3);
    }

    
    std::tuple<unsigned short, int, float, double, uint16_t, uint32_t> pack1;
    std::tuple<Plus, Minus> pack2;
    std::tuple<Square, Identity> pack3;
    //template void iteratePacks(decltype(pack1), decltype(pack2), decltype(pack3));
}

#define TEMPLATE_INSTANTIATION(T, BINARY_OP, UNARY_OP) \
    template T primitives::reduce<T, BINARY_OP>(sycl::queue&, const sycl::span<T>, BINARY_OP); \
    template T primitives::reduce<T, BINARY_OP>(sycl::queue&, const SyclVector<T>&, BINARY_OP); \
    template void primitives::transform<T, UNARY_OP>(sycl::queue&, SyclVector<T>&, UNARY_OP); \
    template void primitives::transform_reduce<T, BINARY_OP, UNARY_OP>(sycl::queue&, sycl::span<T>, BINARY_OP, UNARY_OP); \
    template void primitives::transform_reduce<T, BINARY_OP, UNARY_OP>(sycl::queue&, SyclVector<T>&, BINARY_OP, UNARY_OP); \
    template void primitives::transform<T, UNARY_OP>(sycl::queue&, sycl::span<T>, UNARY_OP); \
    template void primitives::exclusive_scan<T, BINARY_OP>(sycl::queue&, sycl::span<T>, BINARY_OP); \
    template void primitives::exclusive_scan<T, BINARY_OP>(sycl::queue&, SyclVector<T>&, BINARY_OP); \
    template void primitives::inclusive_scan<T, BINARY_OP>(sycl::queue&, sycl::span<T>, BINARY_OP); \
    template void primitives::inclusive_scan<T, BINARY_OP>(sycl::queue&, SyclVector<T>&, BINARY_OP); \
    template void primitives::transform_exclusive_scan<T, BINARY_OP, UNARY_OP>(sycl::queue&, T*, T*, T*, T, BINARY_OP, UNARY_OP);


//TEMPLATE_INSTANTIATION(uint16_t, std::plus<uint16_t>, ConditionFunctor)

template void primitives::exclusive_scan(SyclQueue&, uint16_t*, uint16_t*, uint16_t*, uint16_t, Plus);
template void primitives::exclusive_scan(SyclQueue&, const SyclVector<uint16_t>&, SyclVector<uint16_t>&, uint16_t, Plus);
template void primitives::exclusive_scan(SyclQueue&, const Span<uint16_t>, const Span<uint16_t>, uint16_t, Plus);

template void primitives::exclusive_scan(SyclQueue&, uint16_t*, uint16_t*, uint16_t*, uint16_t, sycl::plus<uint16_t>);
template void primitives::exclusive_scan(SyclQueue&, const SyclVector<uint16_t>&, SyclVector<uint16_t>&, uint16_t, sycl::plus<uint16_t>);
template void primitives::exclusive_scan(SyclQueue&, const Span<uint16_t>, const Span<uint16_t>, uint16_t, sycl::plus<uint16_t>);
