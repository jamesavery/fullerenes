#pragma once
#include <tuple>
#include <array>
#include <fullerenes/sycl-headers/sycl-span.hh>


namespace fullerene_detail {
   
    template<typename... Args1, typename... Args2, std::size_t ... Is>
    auto merge_tuples(std::tuple<Args1&...> t1, std::tuple<Args2&...> t2,
                      std::index_sequence<Is...>) {
        return std::make_tuple(std::tie(std::get<Is>(t1), std::get<Is>(t2))...);
    }

    template<typename Tuple, std::size_t... Is>
    auto rm_const_helper(Tuple&& t, std::index_sequence<Is...>) {
        return std::forward_as_tuple(const_cast<std::remove_const_t<std::tuple_element_t<Is, std::remove_reference_t<Tuple>>>&>(std::get<Is>(std::forward<Tuple>(t)))...);
    }

    template<typename... Args>
    auto rm_const(std::tuple<Args...>&& t) {
        return rm_const_helper(std::move(t), std::index_sequence_for<Args...>{});
    }

    template <std::size_t ...I, typename F> auto with_sequence_impl(F &&func, std::index_sequence<I...>)
    {
        return func(std::integral_constant<std::size_t, I>{}...);
    }

    template <std::size_t N, typename F> auto with_sequence(F &&func)
    {
        return with_sequence_impl(std::forward<F>(func), std::make_index_sequence<N>{});
    }
    
    template<typename... Args1, typename... Args2, std::size_t ... Is, size_t N>
    auto construct_spans_impl(int idx, int count, std::array<int, N>&& size_factors, std::tuple<Args1&...> dst, std::tuple<Args2&...> src,
                        std::index_sequence<Is...>) {
        auto assign_span = [](auto& lhs, auto rhs_ptr, size_t size) {
            lhs = Span(rhs_ptr, size);
        };

        (..., assign_span(std::get<Is>(dst),std::get<Is>(src).data() + size_factors[Is] * ((int)idx), 
                          size_factors[Is]*count   )  );
    }

    template<size_t N, typename... Args1, typename... Args2>
    void construct_spans(int idx, int count, std::array<int, N>&& size_factors, std::tuple<Args1&...> dst, std::tuple<Args2&...> src) {
        construct_spans_impl(idx, count, std::move(size_factors), dst, src, std::make_index_sequence<sizeof...(Args1)>{});
    }
}

template<typename... Args1, typename... Args2>
auto forward_merge_tuples(const std::tuple<Args1&...>& t1, const std::tuple<Args2&...>& t2) {
    static_assert(sizeof...(Args1) == sizeof...(Args2));
    return fullerene_detail::merge_tuples(t1, t2, std::make_index_sequence<sizeof...(Args1)>());
}