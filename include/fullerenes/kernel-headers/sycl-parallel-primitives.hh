#pragma once 
#include <cmath>
#include <numeric>
struct SyclQueue; //Forward declaration
//Unary Operators
struct Identity{
    template <typename T>
    constexpr inline T&& operator()(T&& x) const { return std::forward<T>(x);}

    template <typename T>
    constexpr inline T operator()(const T& x) const { return x;}
};

struct NoOp{
    template <typename T>
    constexpr inline void operator()(const T& x) const {}
};

struct Square{
    template <typename T>
    constexpr inline T operator()(const T& x) const { return x * x;}
};

struct Cube{
    template <typename T>
    constexpr inline T operator()(const T& x) const { return x * x * x;}
};

struct Negate{
    template <typename T>
    constexpr inline T operator()(const T& x) const { return -x;}
};

struct Abs{
    template <typename T>
    constexpr inline T operator()(const T& x) const { return x < 0 ? -x : x;}
};

struct Sign{
    template <typename T>
    constexpr inline T operator()(const T& x) const { return x < 0 ? -1 : (x > 0 ? 1 : 0);}
};

struct Inverse{
    template <typename T>
    constexpr inline T operator()(const T& x) const { return 1 / x;}
};

struct Sqrt{
    template <typename T>
    constexpr inline T operator()(const T& x) const { return std::sqrt(x);}
};

struct Cbrt{
    template <typename T>
    constexpr inline T operator()(const T& x) const { return std::cbrt(x);}
};

struct Exp{
    template <typename T>
    constexpr inline T operator()(const T& x) const { return std::exp(x);}
};

struct Log{
    template <typename T>
    constexpr inline T operator()(const T& x) const { return std::log(x);}
};

struct Log2{
    template <typename T>
    constexpr inline T operator()(const T& x) const { return std::log2(x);}
};

struct Log10{
    template <typename T>
    constexpr inline T operator()(const T& x) const { return std::log10(x);}
};

struct Sin{
    template <typename T>
    constexpr inline T operator()(const T& x) const { return std::sin(x);}
};

struct Cos{
    template <typename T>
    constexpr inline T operator()(const T& x) const { return std::cos(x);}
};

struct Tan{
    template <typename T>
    constexpr inline T operator()(const T& x) const { return std::tan(x);}
};

struct ArcSin{
    template <typename T>
    constexpr inline T operator()(const T& x) const { return std::asin(x);}
};

struct ArcCos{
    template <typename T>
    constexpr inline T operator()(const T& x) const { return std::acos(x);}
};

struct ArcTan{
    template <typename T>
    constexpr inline T operator()(const T& x) const { return std::atan(x);}
};

struct Sinh{
    template <typename T>
    constexpr inline T operator()(const T& x) const { return std::sinh(x);}
};

struct Cosh{
    template <typename T>
    constexpr inline T operator()(const T& x) const { return std::cosh(x);}
};

struct Tanh{
    template <typename T>
    constexpr inline T operator()(const T& x) const { return std::tanh(x);}
};

struct ArcSinh{
    template <typename T>
    constexpr inline T operator()(const T& x) const { return std::asinh(x);}
};

struct ArcCosh{
    template <typename T>
    constexpr inline T operator()(const T& x) const { return std::acosh(x);}
};

struct ArcTanh{
    template <typename T>
    constexpr inline T operator()(const T& x) const { return std::atanh(x);}
};

struct Erf{
    template <typename T>
    constexpr inline T operator()(const T& x) const { return std::erf(x);}
};

struct Erfc{
    template <typename T>
    constexpr inline T operator()(const T& x) const { return std::erfc(x);}
};

struct Tgamma{
    template <typename T>
    constexpr inline T operator()(const T& x) const { return std::tgamma(x);}
};

struct Lgamma{
    template <typename T>
    constexpr inline T operator()(const T& x) const { return std::lgamma(x);}
};

struct Round{
    template <typename T>
    constexpr inline T operator()(const T& x) const { return std::round(x);}
};

struct AbsDiff{
    template <typename T>
    constexpr inline T operator()(const T& x, const T& y) const { return x > y ? x - y : y - x;}
};

struct Hypot{
    template <typename T>
    constexpr inline T operator()(const T& x, const T& y) const { return std::hypot(x, y);}
};

struct Floor{
    template <typename T>
    constexpr inline T operator()(const T& x) const { return std::floor(x);}
};

struct Ceil{
    template <typename T>
    constexpr inline T operator()(const T& x) const { return std::ceil(x);}
};




//Associative Binary Operators
struct Plus{
    template <typename T>
    constexpr inline T operator()(const T& x, const T& y) const { return x + y;}

    template <typename T, typename U, typename = std::enable_if_t<!std::is_same_v<T, U>>>
    constexpr inline T operator()(const T& x, const U& y) const { return x + y;}
};

struct Multiply{
    template <typename T>
    constexpr inline T operator()(const T& x, const T& y) const { return x * y;}

    template <typename T, typename U, typename = std::enable_if_t<!std::is_same_v<T, U>>>
    constexpr inline T operator()(const T& x, const U& y) const { return x * y;}
};

struct Max{
    template <typename T>
    constexpr inline T operator()(const T& x, const T& y) const { return x > y ? x : y;}

    template <typename T, typename U, typename = std::enable_if_t<!std::is_same_v<T, U>>>
    constexpr inline T operator()(const T& x, const U& y) const { return x > y ? x : y;}
};

struct Min{
    template <typename T>
    constexpr inline T operator()(const T& x, const T& y) const { return x < y ? x : y;}

    template <typename T, typename U, typename = std::enable_if_t<!std::is_same_v<T, U>>>
    constexpr inline T operator()(const T& x, const U& y) const { return x < y ? x : y;}
};

struct GreatestCommonDivisor{
    template <typename T>
    constexpr inline T operator()(const T& x, const T& y) const { return std::gcd(x, y);}

    template <typename T, typename U>
    constexpr inline T operator()(const T& x, const U& y) const { return std::gcd(x, y);}

    template <typename T, typename U>
    constexpr inline T operator()(const U& x, const T& y) const { return std::gcd(x, y);}
};

struct LeastCommonMultiple{
    template <typename T>
    constexpr inline T operator()(const T& x, const T& y) const { return std::lcm(x, y);}

    template <typename T, typename U>
    constexpr inline T operator()(const T& x, const U& y) const { return std::lcm(x, y);}

    template <typename T, typename U>
    constexpr inline T operator()(const U& x, const T& y) const { return std::lcm(x, y);}
};

//Unary Predicates

struct AlwaysTrue{
    template <typename T>
    constexpr inline bool operator()(const T& x) const { return true;}
};

struct AlwaysFalse{
    template <typename T>
    constexpr inline bool operator()(const T& x) const { return false;}
};

struct IsTrue{
    template <typename T>
    constexpr inline bool operator()(const T& x) const { return x;}
};

struct IsFalse{
    template <typename T>
    constexpr inline bool operator()(const T& x) const { return !x;}
};

template <typename T>
struct IsEqualTo{
    const T value;
    constexpr IsEqualTo(const T& value) : value(value) {}
    constexpr inline bool operator()(const T& x) const { return x == value;}
};

template <typename T>
struct IsNotEqualTo{
    const T value;
    constexpr IsNotEqualTo(const T& value) : value(value) {}
    constexpr inline bool operator()(const T& x) const { return x != value;}
};

template <typename T>
struct IsLessThan{
    const T value;
    constexpr IsLessThan(const T& value) : value(value) {}
    constexpr inline bool operator()(const T& x) const { return x < value;}
};

template <typename T>
struct IsGreaterThan{
    const T value;
    constexpr IsGreaterThan(const T& value) : value(value) {}
    constexpr inline bool operator()(const T& x) const { return x > value;}
};

struct IsPositive{
    template <typename T>
    constexpr inline bool operator()(const T& x) const { return x > 0;}
};

struct IsNegative{
    template <typename T>
    constexpr inline bool operator()(const T& x) const { return x < 0;}
};

struct IsOdd{
    template <typename T>
    constexpr inline bool operator()(const T& x) const { return x % 2;}
};

struct IsEven{
    template <typename T>
    constexpr inline bool operator()(const T& x) const { return !(x % 2);}
};

struct IsFinite{
    template <typename T>
    constexpr inline bool operator()(const T& x) const { return std::isfinite(x);}
};

struct IsInf{
    template <typename T>
    constexpr inline bool operator()(const T& x) const { return std::isinf(x);}
};

struct IsNan{
    template <typename T>
    constexpr inline bool operator()(const T& x) const { return std::isnan(x);}
};

struct IsNormal{
    template <typename T>
    constexpr inline bool operator()(const T& x) const { return std::isnormal(x);}
};

struct IsZero{
    template <typename T>
    constexpr inline bool operator()(const T& x) const { return x == 0;}
};

struct IsNonZero{
    template <typename T>
    constexpr inline bool operator()(const T& x) const { return x != 0;}
};

struct IsNonNegative{
    template <typename T>
    constexpr inline bool operator()(const T& x) const { return x >= 0;}
};

struct IsNonPositive{
    template <typename T>
    constexpr inline bool operator()(const T& x) const { return x <= 0;}
};



//Binary Predicates
struct Equal{
    template <typename T>
    constexpr inline bool operator()(const T& x, const T& y) const { return x == y;}
};

struct NotEqual{
    template <typename T>
    constexpr inline bool operator()(const T& x, const T& y) const { return x != y;}
};

struct Less{
    template <typename T>
    constexpr inline bool operator()(const T& x, const T& y) const { return x < y;}
};

struct Greater{
    template <typename T>
    constexpr inline bool operator()(const T& x, const T& y) const { return x > y;}
};

struct LessEqual{
    template <typename T>
    constexpr inline bool operator()(const T& x, const T& y) const { return x <= y;}
};

struct GreaterEqual{
    template <typename T>
    constexpr inline bool operator()(const T& x, const T& y) const { return x >= y;}
};

struct LogicalAnd{
    template <typename T>
    constexpr inline bool operator()(const T& x, const T& y) const { return x && y;}
};

struct LogicalOr{
    template <typename T>
    constexpr inline bool operator()(const T& x, const T& y) const { return x || y;}
};

struct LogicalNot{
    template <typename T>
    constexpr inline bool operator()(const T& x) const { return !x;}
};

struct FloatEqual{
    const float eps = std::numeric_limits<float>::epsilon() * 1e2;
    FloatEqual() = default;
    FloatEqual(float eps) : eps(eps) {}

    template <typename T>
    constexpr inline bool operator()(const T& x, const T& y) const {
        T diff = std::abs(x - y);
        T max_v = std::max(std::abs(x), std::abs(y));
        return (diff / (max_v > eps ? max_v : 1)) < eps;
    }
};

namespace primitives {
    // Segmented algorithms retained during std::algorithm migration.
    template <typename InputKeys, typename InputValues, typename OutputKeys, typename OutputValues, typename BinaryPredicate = Equal, typename BinaryOp = Plus>
    inline void reduce_by_segment(SyclQueue& Q, InputKeys&& keys, InputValues&& values, OutputKeys&& out_keys, OutputValues&& out_values, BinaryPredicate pred = BinaryPredicate{}, BinaryOp op = BinaryOp{}) {
        (void)Q;
        auto k_it = keys.begin();
        auto v_it = values.begin();
        auto k_end = keys.end();
        auto ok_it = out_keys.begin();
        auto ov_it = out_values.begin();

        if (k_it == k_end) {
            return;
        }

        auto current_key = *k_it;
        auto accumulated = *v_it;
        ++k_it;
        ++v_it;

        while (k_it != k_end) {
            if (!pred(*k_it, current_key)) {
                *ok_it++ = current_key;
                *ov_it++ = accumulated;
                current_key = *k_it;
                accumulated = *v_it;
            } else {
                accumulated = op(accumulated, *v_it);
            }
            ++k_it;
            ++v_it;
        }

        *ok_it = current_key;
        *ov_it = accumulated;
    }
}