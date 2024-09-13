
#if __INTEL_LLVM_COMPILER > 20240000
    #define ONEDPL_USE_PREDEFINED_POLICIES 0
    #include <oneapi/dpl/algorithm>
    #include <oneapi/dpl/execution>
    #include <oneapi/dpl/iterator>
    #include <oneapi/dpl/numeric>
#endif

#include <fullerenes/sycl-headers/sycl-span.hh>
#include <fullerenes/sycl-headers/sycl-vector.hh>
#include <fullerenes/sycl-headers/sycl-device-queue.hh>
#include <fullerenes/sycl-headers/sycl-parallel-primitives.hh>
#include <fullerenes/sycl-headers/sycl-status-enum.hh>
#include <execution>
#include "queue-impl.cc"    

namespace primitives{


#if __INTEL_LLVM_COMPILER > 20240000
    #define BEGIN(x) (static_cast< std::decay_t<decltype(x[0])>* >(x.begin()))
    #define END(x) (static_cast< std::decay_t<decltype(x[0])>* >(x.end()))

    template <typename InputContainer, typename OutputContainer, typename Init, typename BinaryOp>
    void inline exclusive_scan(SyclQueue& Q, InputContainer&& in_vec, OutputContainer&& out_vec, Init init, BinaryOp op) {
        auto policy = oneapi::dpl::execution::make_device_policy((*Q));
        oneapi::dpl::exclusive_scan(policy, BEGIN(in_vec), END(in_vec), BEGIN(out_vec), init, op);
    }

    template <typename InputContainer, typename OutputContainer, typename Init, typename BinaryOp>
    void inline inclusive_scan(SyclQueue& Q, InputContainer&& in_vec, OutputContainer&& out_vec, Init init, BinaryOp op) {
        auto policy = oneapi::dpl::execution::make_device_policy((*Q));
        oneapi::dpl::inclusive_scan(policy, BEGIN(in_vec), END(in_vec), BEGIN(out_vec), op, init);
    }

    template <typename InputContainer, typename Init, typename BinaryOp>
    auto inline reduce(SyclQueue& Q, InputContainer&& vec, Init init, BinaryOp op) -> decltype(init) {
        auto policy = oneapi::dpl::execution::make_device_policy((*Q));
        return oneapi::dpl::reduce(policy, BEGIN(vec), END(vec), init, op);
    }

    template <typename InputContainer, typename OutputContainer, typename Init, typename BinaryOp, typename UnaryOp>
    void inline transform_inclusive_scan(SyclQueue& Q, InputContainer&& in_vec, OutputContainer&& out_vec, Init init, BinaryOp binary_op, UnaryOp unary_op) {
        auto policy = oneapi::dpl::execution::make_device_policy((*Q));
        oneapi::dpl::transform_inclusive_scan(policy, BEGIN(in_vec), END(in_vec), BEGIN(out_vec), binary_op, unary_op, init);
    }

    template <typename InputContainer1, typename InputContainer2, typename OutputContainer, typename Init, typename ReduceOp, typename TransformOp>
    void inline transform_inclusive_scan(SyclQueue& Q, InputContainer1&& vec1, InputContainer2&& vec2, OutputContainer&& out_vec, Init init, ReduceOp reduce_op, TransformOp transform_op) {
        auto policy = oneapi::dpl::execution::make_device_policy((*Q));
        oneapi::dpl::transform_inclusive_scan(policy, BEGIN(vec1), END(vec1), BEGIN(vec2), BEGIN(out_vec), reduce_op, transform_op, init);
    }

    template <typename InputContainer, typename OutputContainer, typename Init, typename BinaryOp, typename UnaryOp>
    void inline transform_exclusive_scan(SyclQueue& Q, InputContainer&& in_vec, OutputContainer&& out_vec, Init init, BinaryOp binary_op, UnaryOp unary_op) {
        auto policy = oneapi::dpl::execution::make_device_policy((*Q));
        oneapi::dpl::transform_exclusive_scan(policy, BEGIN(in_vec), END(in_vec), BEGIN(out_vec), init, binary_op, unary_op);
    }

    template <typename InputContainer1, typename InputContainer2, typename OutputContainer, typename Init, typename ReduceOp, typename TransformOp>
    void inline transform_exclusive_scan(SyclQueue& Q, InputContainer1&& vec1, InputContainer2&& vec2, OutputContainer&& out_vec, Init init, ReduceOp reduce_op, TransformOp transform_op) {
        auto policy = oneapi::dpl::execution::make_device_policy((*Q));
        oneapi::dpl::transform_exclusive_scan(policy, BEGIN(vec1), END(vec1), BEGIN(vec2), BEGIN(out_vec), init, reduce_op, transform_op);
    }

    template <typename InputContainer, typename OutputContainer, typename UnaryOp>
    void inline transform(SyclQueue& Q, InputContainer&& in_vec, OutputContainer&& out_vec, UnaryOp unary_op) {
        auto policy = oneapi::dpl::execution::make_device_policy((*Q));
        oneapi::dpl::transform(policy, BEGIN(in_vec), END(in_vec), BEGIN(out_vec), unary_op);
    }

    template <typename InputContainer1, typename InputContainer2, typename OutputContainer, typename BinaryOp>
    void inline transform(SyclQueue& Q, InputContainer1&& vec1, InputContainer2&& vec2, OutputContainer&& out_vec, BinaryOp binary_op) {
        auto policy = oneapi::dpl::execution::make_device_policy((*Q));
        oneapi::dpl::transform(policy, BEGIN(vec1), END(vec1), BEGIN(vec2), BEGIN(out_vec), binary_op);
    }

    template <typename InputContainer, typename Init, typename BinaryOp, typename UnaryOp>
    auto inline transform_reduce(SyclQueue& Q, InputContainer&& vec, Init init, BinaryOp binary_op, UnaryOp unary_op) -> decltype(init) {
        auto policy = oneapi::dpl::execution::make_device_policy((*Q));
        return oneapi::dpl::transform_reduce(policy, BEGIN(vec), END(vec), init, binary_op, unary_op);
    }

    template <typename InputContainer1, typename InputContainer2, typename Init, typename ReduceOp, typename TransformOp>
    auto inline transform_reduce(SyclQueue& Q, InputContainer1&& vec1, InputContainer2&& vec2, Init init, ReduceOp reduce_op, TransformOp transform_op) -> decltype(init) {
        auto policy = oneapi::dpl::execution::make_device_policy((*Q));
        return oneapi::dpl::transform_reduce(policy, BEGIN(vec1), END(vec1), BEGIN(vec2), init, reduce_op, transform_op);
    }

    //any_of
    template <typename InputContainer, typename Predicate>
    bool inline any_of(SyclQueue& Q, InputContainer&& vec, Predicate pred) {
        auto policy = oneapi::dpl::execution::make_device_policy((*Q));
        return oneapi::dpl::any_of(policy, BEGIN(vec), END(vec), pred);
    }

    //all_of
    template <typename InputContainer, typename Predicate>
    bool inline all_of(SyclQueue& Q, InputContainer&& vec, Predicate pred) {
        auto policy = oneapi::dpl::execution::make_device_policy((*Q));
        return oneapi::dpl::all_of(policy, BEGIN(vec), END(vec), pred);
    }

    //none_of
    template <typename InputContainer, typename Predicate>
    bool inline none_of(SyclQueue& Q, InputContainer&& vec, Predicate pred) {
        auto policy = oneapi::dpl::execution::make_device_policy((*Q));
        return oneapi::dpl::none_of(policy, BEGIN(vec), END(vec), pred);
    }

    //for_each
    template <typename InputContainer, typename Function>
    void inline for_each(SyclQueue& Q, InputContainer&& vec, Function func) {
        auto policy = oneapi::dpl::execution::make_device_policy((*Q));
        oneapi::dpl::for_each(policy, BEGIN(vec), END(vec), func);
    }

    //for_each_n
    template <typename InputContainer, typename Function>
    void inline for_each_n(SyclQueue& Q, InputContainer&& vec, size_t n, Function func) {
        auto policy = oneapi::dpl::execution::make_device_policy((*Q));
        oneapi::dpl::for_each_n(policy, BEGIN(vec), n, func);
    }

    //count
    template <typename InputContainer, typename T>
    size_t inline count(SyclQueue& Q, InputContainer&& vec, T value) {
        auto policy = oneapi::dpl::execution::make_device_policy((*Q));
        return oneapi::dpl::count(policy, BEGIN(vec), END(vec), value);
    }

    //count_if
    template <typename InputContainer, typename Predicate>
    size_t inline count_if(SyclQueue& Q, InputContainer&& vec, Predicate pred) {
        auto policy = oneapi::dpl::execution::make_device_policy((*Q));
        return oneapi::dpl::count_if(policy, BEGIN(vec), END(vec), pred);
    }

    //copy
    template <typename InputContainer, typename OutputContainer>
    void inline copy(SyclQueue& Q, InputContainer&& in_vec, OutputContainer&& out_vec) {
        auto policy = oneapi::dpl::execution::make_device_policy((*Q));
        oneapi::dpl::copy(policy, BEGIN(in_vec), END(in_vec), BEGIN(out_vec));
    }

    //copy_if
    template <typename InputContainer, typename OutputContainer, typename Predicate>
    void inline copy_if(SyclQueue& Q, InputContainer&& in_vec, OutputContainer&& out_vec, Predicate unary_op) {
        auto policy = oneapi::dpl::execution::make_device_policy((*Q));
        oneapi::dpl::copy_if(policy, BEGIN(in_vec), END(in_vec), BEGIN(out_vec), unary_op);
    }

    //copy_n
    template <typename InputContainer, typename OutputContainer>
    void inline copy_n(SyclQueue& Q, InputContainer&& in_vec, size_t n, OutputContainer&& out_vec) {
        auto policy = oneapi::dpl::execution::make_device_policy((*Q));
        oneapi::dpl::copy_n(policy, BEGIN(in_vec), n, BEGIN(out_vec));
    }

    //fill
    template <typename InputContainer, typename T>
    void inline fill(SyclQueue& Q, InputContainer&& vec, T value) {
        auto policy = oneapi::dpl::execution::make_device_policy((*Q));
        oneapi::dpl::fill(policy, BEGIN(vec), END(vec), value);
    }

    //sort
    template <typename InputContainer, typename BinaryPredicate>
    void inline sort(SyclQueue& Q, InputContainer&& vec, BinaryPredicate binary_op) {
        auto policy = oneapi::dpl::execution::make_device_policy((*Q));
        oneapi::dpl::sort(policy, BEGIN(vec), END(vec), binary_op);
    }

    //iota
    template <typename InputContainer, typename T>
    void inline iota(SyclQueue& Q, InputContainer&& vec, T value) {
        auto policy = oneapi::dpl::execution::make_device_policy((*Q));
        auto begin = BEGIN(vec);
        oneapi::dpl::for_each(policy, BEGIN(vec), END(vec), [begin, value](auto& x) { auto index = &x - begin; x = value + index; });
    }

    //OneAPI Algorithms

    //inclusive_scan_by_segment
    template <typename InputKeys, typename InputValues, typename OutputValues, typename BinaryPredicate, typename BinaryOp>
    void inline inclusive_scan_by_segment(SyclQueue& Q, InputKeys&& keys, InputValues&& values, OutputValues&& out_values, BinaryPredicate pred, BinaryOp op) {
        auto policy = oneapi::dpl::execution::make_device_policy((*Q));
        oneapi::dpl::inclusive_scan_by_segment(policy, BEGIN(keys), END(keys), BEGIN(values), BEGIN(out_values), pred, op);
    }

    //exclusive_scan_by_segment
    template <typename InputKeys, typename InputValues, typename OutputValues, typename BinaryPredicate, typename BinaryOp>
    void inline exclusive_scan_by_segment(SyclQueue& Q, InputKeys&& keys, InputValues&& values, OutputValues&& out_values, BinaryPredicate pred, BinaryOp op) {
        auto policy = oneapi::dpl::execution::make_device_policy((*Q));
        oneapi::dpl::exclusive_scan_by_segment(policy, BEGIN(keys), END(keys), BEGIN(values), BEGIN(out_values), pred, op);
    }

    //reduce_by_segment
    template <typename InputKeys, typename InputValues, typename OutputKeys, typename OutputValues, typename BinaryPredicate, typename BinaryOp>
    void inline reduce_by_segment(SyclQueue& Q, InputKeys&& keys, InputValues&& values, OutputKeys&& out_keys, OutputValues&& out_values, BinaryPredicate pred, BinaryOp op) {
        auto policy = oneapi::dpl::execution::make_device_policy((*Q));
        oneapi::dpl::reduce_by_segment(policy, BEGIN(keys), END(keys), BEGIN(values), BEGIN(out_keys), BEGIN(out_values), pred, op);
    }

    //histogram


#else

    //Lets define a SFINAE helper to check if the type can be called with T1 and T2
    #define BEGIN(x) (static_cast< std::decay_t<decltype(x[0])>* >(x.begin()))
    #define END(x) (static_cast< std::decay_t<decltype(x[0])>* >(x.end()))

    template <typename InputContainer, typename OutputContainer, typename Init, typename BinaryOp>
    void inline exclusive_scan(SyclQueue& Q, InputContainer&& in_vec, OutputContainer&& out_vec, Init init, BinaryOp op) {
        Q.wait();
        std::exclusive_scan(std::execution::par_unseq, BEGIN(in_vec), END(in_vec), BEGIN(out_vec), init, op);
    }

    template <typename InputContainer, typename OutputContainer, typename Init, typename BinaryOp>
    void inline inclusive_scan(SyclQueue& Q, InputContainer&& in_vec, OutputContainer&& out_vec, Init init, BinaryOp op) {
        Q.wait();
        std::inclusive_scan(std::execution::par_unseq, BEGIN(in_vec), END(in_vec), BEGIN(out_vec), op, init);
    }

    template <typename InputContainer, typename Init, typename BinaryOp>
    auto inline reduce(SyclQueue& Q, InputContainer&& vec, Init init, BinaryOp op) -> decltype(init) {
        Q.wait();
        return std::reduce(std::execution::par_unseq, BEGIN(vec), END(vec), init, op);
    }

    template <typename InputContainer, typename OutputContainer, typename Init, typename BinaryOp, typename UnaryOp>
    void inline transform_inclusive_scan(SyclQueue& Q, InputContainer&& in_vec, OutputContainer&& out_vec, Init init, BinaryOp binary_op, UnaryOp unary_op) {
        Q.wait();
        std::transform_inclusive_scan(std::execution::par_unseq, BEGIN(in_vec), END(in_vec), BEGIN(out_vec), binary_op, unary_op, init);
    }

    template <typename InputContainer1, typename InputContainer2, typename OutputContainer, typename Init, typename ReduceOp, typename TransformOp>
    void inline transform_inclusive_scan(SyclQueue& Q, InputContainer1&& vec1, InputContainer2&& vec2, OutputContainer&& out_vec, Init init, ReduceOp reduce_op, TransformOp transform_op) {
        Q.wait();
        std::transform_inclusive_scan(std::execution::par_unseq, BEGIN(vec1), END(vec1), BEGIN(vec2), BEGIN(out_vec), reduce_op, transform_op, init);
    }

    template <typename InputContainer, typename OutputContainer, typename Init, typename BinaryOp, typename UnaryOp>
    void inline transform_exclusive_scan(SyclQueue& Q, InputContainer&& in_vec, OutputContainer&& out_vec, Init init, BinaryOp binary_op, UnaryOp unary_op) {
        Q.wait();
        std::transform_exclusive_scan(std::execution::par_unseq, BEGIN(in_vec), END(in_vec), BEGIN(out_vec), init, binary_op, unary_op);
    }

    template <typename InputContainer1, typename InputContainer2, typename OutputContainer, typename Init, typename ReduceOp, typename TransformOp>
    void inline transform_exclusive_scan(SyclQueue& Q, InputContainer1&& vec1, InputContainer2&& vec2, OutputContainer&& out_vec, Init init, ReduceOp reduce_op, TransformOp transform_op) {
        Q.wait();
        std::transform_exclusive_scan(std::execution::par_unseq, BEGIN(vec1), END(vec1), BEGIN(vec2), BEGIN(out_vec), init, reduce_op, transform_op);
    }

    template <typename InputContainer, typename OutputContainer, typename UnaryOp>
    void inline transform(SyclQueue& Q, InputContainer&& in_vec, OutputContainer&& out_vec, UnaryOp unary_op) {
        Q.wait();
        std::transform(std::execution::par_unseq, BEGIN(in_vec), END(in_vec), BEGIN(out_vec), unary_op);
    }

    template <typename InputContainer1, typename InputContainer2, typename OutputContainer, typename BinaryOp>
    void inline transform(SyclQueue& Q, InputContainer1&& vec1, InputContainer2&& vec2, OutputContainer&& out_vec, BinaryOp binary_op) {
        Q.wait();
        std::transform(std::execution::par_unseq, BEGIN(vec1), END(vec1), BEGIN(vec2), BEGIN(out_vec), binary_op);
    }

    template <typename InputContainer, typename Init, typename BinaryOp, typename UnaryOp>
    auto inline transform_reduce(SyclQueue& Q, InputContainer&& vec, Init init, BinaryOp binary_op, UnaryOp unary_op) -> decltype(init) {
        Q.wait();
        return std::transform_reduce(std::execution::par_unseq, BEGIN(vec), END(vec), init, binary_op, unary_op);
    }

    template <typename InputContainer1, typename InputContainer2, typename Init, typename ReduceOp, typename TransformOp>
    auto inline transform_reduce(SyclQueue& Q, InputContainer1&& vec1, InputContainer2&& vec2, Init init, ReduceOp reduce_op, TransformOp transform_op) -> decltype(init) {
        Q.wait();
        return std::transform_reduce(std::execution::par_unseq, BEGIN(vec1), END(vec1), BEGIN(vec2), init, reduce_op, transform_op);
    }

    //any_of
    template <typename InputContainer, typename Predicate>
    bool inline any_of(SyclQueue& Q, InputContainer&& vec, Predicate pred) {
        Q.wait();
        return std::any_of(std::execution::par_unseq, BEGIN(vec), END(vec), pred);
    }

    //all_of
    template <typename InputContainer, typename Predicate>
    bool inline all_of(SyclQueue& Q, InputContainer&& vec, Predicate pred) {
        Q.wait();
        return std::all_of(std::execution::par_unseq, BEGIN(vec), END(vec), pred);
    }

    //none_of
    template <typename InputContainer, typename Predicate>
    bool inline none_of(SyclQueue& Q, InputContainer&& vec, Predicate pred) {
        Q.wait();
        return std::none_of(std::execution::par_unseq, BEGIN(vec), END(vec), pred);
    }

    //for_each
    template <typename InputContainer, typename Function>
    void inline for_each(SyclQueue& Q, InputContainer&& vec, Function func) {
        Q.wait();
        std::for_each(std::execution::par_unseq, BEGIN(vec), END(vec), func);
    }

    //for_each_n
    template <typename InputContainer, typename Function>
    void inline for_each_n(SyclQueue& Q, InputContainer&& vec, size_t n, Function func) {
        Q.wait();
        std::for_each_n(std::execution::par_unseq, BEGIN(vec), n, func);
    }

    //count
    template <typename InputContainer, typename T>
    size_t inline count(SyclQueue& Q, InputContainer&& vec, T value) {
        Q.wait();
        return std::count(std::execution::par_unseq, BEGIN(vec), END(vec), value);
    }

    //count_if
    template <typename InputContainer, typename Predicate>
    size_t inline count_if(SyclQueue& Q, InputContainer&& vec, Predicate pred) {
        Q.wait();
        return std::count_if(std::execution::par_unseq, BEGIN(vec), END(vec), pred);
    }

    //copy
    template <typename InputContainer, typename OutputContainer>
    void inline copy(SyclQueue& Q, InputContainer&& in_vec, OutputContainer&& out_vec) {
        Q.wait();
        std::copy(std::execution::par_unseq, BEGIN(in_vec), END(in_vec), BEGIN(out_vec));
    }

    //copy_if
    template <typename InputContainer, typename OutputContainer, typename Predicate>
    void inline copy_if(SyclQueue& Q, InputContainer&& in_vec, OutputContainer&& out_vec, Predicate unary_op) {
        Q.wait();
        std::copy_if(std::execution::par_unseq, BEGIN(in_vec), END(in_vec), BEGIN(out_vec), unary_op);
    }

    //copy_n
    template <typename InputContainer, typename OutputContainer>
    void inline copy_n(SyclQueue& Q, InputContainer&& in_vec, size_t n, OutputContainer&& out_vec) {
        Q.wait();
        std::copy_n(std::execution::par_unseq, BEGIN(in_vec), n, BEGIN(out_vec));
    }

    //fill
    template <typename InputContainer, typename T>
    void inline fill(SyclQueue& Q, InputContainer&& vec, T value) {
        Q.wait();
        std::fill(std::execution::par_unseq, BEGIN(vec), END(vec), value);
    }

    //sort
    template <typename InputContainer, typename BinaryPredicate>
    void inline sort(SyclQueue& Q, InputContainer&& vec, BinaryPredicate binary_op) {
        Q.wait();
        std::sort(std::execution::par_unseq, BEGIN(vec), END(vec), binary_op);
    }
#endif

}

