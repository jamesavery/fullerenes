#include "primitives.cc"
namespace primitives{

    template <typename T1, typename T2, template <typename> class ContainerType, typename BinaryOp, typename UnaryOp, typename BinaryPredicate, typename Predicate>
    void myFunction() {
        SyclQueue Q;
        ContainerType<T1> input;
        ContainerType<T2> output;
        T1 res;
        res = reduce(Q, input, T1{}, BinaryOp{});
        res = transform_reduce(Q, input, T1{}, BinaryOp{}, UnaryOp{});

        exclusive_scan(Q, input, output, T1{}, BinaryOp{});
        inclusive_scan(Q, input, output, T1{}, BinaryOp{});
        transform_exclusive_scan(Q, input, output, T1{}, BinaryOp{}, UnaryOp{});
        transform_inclusive_scan(Q, input, output, T1{}, BinaryOp{}, UnaryOp{});
        transform(Q, input, output, UnaryOp{});
        copy_if(Q, input, output, Predicate{});
        fill(Q, input, T1{});
        sort(Q, input, BinaryPredicate{});
        any_of(Q, input, Predicate{});
        all_of(Q, input, Predicate{});
        none_of(Q, input, Predicate{});
        fill(Q, input, T1{});
    }

    /* 
    using Types = std::tuple<uint16_t, float>;
    using UnaryOperators = std::tuple<Identity, Negate, Square, Cube>; 
    using BinaryOperators = std::tuple<Plus, Multiply, Min, Max>;
    using BinaryPredicates = std::tuple<Less, Greater, Equal, NotEqual, GreaterEqual, LessEqual>;
    using Predicates = std::tuple<IsTrue, IsFalse>;
    */
    using Types = std::tuple<float>;
    using UnaryOperators = std::tuple<Identity>; 
    using BinaryOperators = std::tuple<Plus>;
    using BinaryPredicates = std::tuple<Equal>;
    using Predicates = std::tuple<IsTrue>;


    template <typename Tuple, std::size_t Index>
    using tuple_element_t = typename std::tuple_element<Index, Tuple>::type;

    template <std::size_t I, std::size_t J, std::size_t K, std::size_t L, std::size_t M, std::size_t N>
    struct InstantiateFunctions {
        static void instantiate() {
            using T1 = tuple_element_t<Types, I>;
            using T2 = tuple_element_t<Types, J>;
            using UnaryOp = tuple_element_t<UnaryOperators, K>;
            using BinaryOp = tuple_element_t<BinaryOperators, L>;
            using BinaryPredicate = tuple_element_t<BinaryPredicates, M>;
            using Predicate = tuple_element_t<Predicates, N>;

            myFunction<T1, T2, Span,        BinaryOp, UnaryOp, BinaryPredicate, Predicate>();
            myFunction<T1, T2, SyclVector,  BinaryOp, UnaryOp, BinaryPredicate, Predicate>();

            if constexpr (N + 1 < std::tuple_size_v<Predicates>){
                InstantiateFunctions<I, J, K, L, M, N + 1>::instantiate();
            }else if constexpr (M + 1 < std::tuple_size_v<BinaryPredicates>) {
                InstantiateFunctions<I, J, K, L, M + 1, 0>::instantiate();
            } else if constexpr (L + 1 < std::tuple_size_v<BinaryOperators>) {
                InstantiateFunctions<I, J, K, L + 1, 0, 0>::instantiate();
            } else if constexpr (K + 1 < std::tuple_size_v<UnaryOperators>) {
                InstantiateFunctions<I, J, K + 1, 0, 0, 0>::instantiate();
            } else if constexpr (J + 1 < std::tuple_size_v<Types>) {
                InstantiateFunctions<I, J + 1, 0, 0, 0, 0>::instantiate();
            } else if constexpr (I + 1 < std::tuple_size_v<Types>) {
                InstantiateFunctions<I + 1, 0, 0, 0, 0, 0>::instantiate();
            }
        }
    };
    



    void instantiateAllFunctions() {
        if constexpr (std::tuple_size_v<Types> == 0) return;
        else {
            InstantiateFunctions<0, 0, 0, 0, 0, 0>::instantiate();
        }
    }
}