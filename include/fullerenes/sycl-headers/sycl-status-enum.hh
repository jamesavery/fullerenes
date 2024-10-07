#pragma once
#include <type_traits>
#include <iostream>
#include <array>

enum class StatusEnum
{
    EMPTY = 1 << 0,            // 0000 0001
    CONVERGED_2D = 1 << 1,     // 0000 0010
    CONVERGED_3D = 1 << 2,     // 0000 0100
    PLZ_CHECK = 1 << 3,        // 0000 1000
    FAILED_2D = 1 << 4,        // 0001 0000
    FAILED_3D = 1 << 5,        // 0010 0000
    DUAL_INITIALIZED = 1 << 6, // 0100 0000
    CUBIC_INITIALIZED = 1 << 7 // 1000 0000
};

struct StatusFlag
{
    using enum StatusEnum;

    StatusFlag() : flag_(0) {}
    StatusFlag(const StatusFlag& flag) = default;
    StatusFlag& operator=(const StatusFlag& flag) = default;
    StatusFlag(StatusFlag&& flag) = default;
    StatusFlag& operator=(StatusFlag&& flag) = default;

    StatusFlag& operator=(StatusEnum flag) { flag_ = static_cast<int>(flag); return *this; }

    template <typename T>
    StatusFlag& operator=(const T& flag) { flag_ = static_cast<int>(flag); return *this; }

    template <typename T>
    StatusFlag(T flag) : flag_(static_cast<int>(flag)) {}

    template <typename T>
    StatusFlag operator|=(T flag) { return flag_ |= static_cast<int>(flag); }
    template <typename T>
    StatusFlag operator&=(T flag) { return flag_ &= static_cast<int>(flag); }
    template <typename T>
    StatusFlag operator^=(T flag) { return flag_ ^= static_cast<int>(flag); }

    template <typename T>
    StatusFlag operator|(T flag) const { return flag_ | static_cast<int>(flag); }
    template <typename T>
    StatusFlag operator&(T flag) const { return flag_ & static_cast<int>(flag); }
    template <typename T>
    StatusFlag operator^(T flag) const { return flag_ ^ static_cast<int>(flag); }





    operator int() const { return static_cast<int>(flag_); }
    operator StatusEnum() const { return static_cast<StatusEnum>(flag_); }

    template <typename T>
    inline constexpr bool operator==(T flag) const { return flag_ == static_cast<int>(flag); }

    template <typename T>
    inline constexpr bool operator!=(T flag) const { return flag_ != static_cast<int>(flag); }

    template <typename T>
    inline constexpr bool operator<(T flag) const { return flag_ < static_cast<int>(flag); }

    template <typename T>
    inline constexpr bool operator>(T flag) const { return flag_ > static_cast<int>(flag); }

    template <typename T>
    inline constexpr bool operator<=(T flag) const { return flag_ <= static_cast<int>(flag); }

    template <typename T>
    inline constexpr bool operator>=(T flag) const { return flag_ >= static_cast<int>(flag); }

    template <typename T>
    inline constexpr bool operator&&(T flag) const { return flag_ && static_cast<int>(flag); }

    template <typename T>
    inline constexpr bool operator||(T flag) const { return flag_ || static_cast<int>(flag); }

    template <typename T>
    inline constexpr bool operator!() const { return !flag_; }

    template <typename T>
    inline constexpr bool is_set(T flag) const { return flag_ & static_cast<int>(flag); }

    template <typename T>
    inline constexpr bool is_not_set(T flag) const { return !(flag_ & static_cast<int>(flag)); }

    template <typename T>
    inline constexpr void set(T flag) { flag_ |= static_cast<int>(flag); }

    template <typename T>
    inline constexpr void unset(T flag) { flag_ &= ~static_cast<int>(flag); }

    template <typename T>
    inline constexpr void toggle(T flag) { flag_ ^= static_cast<int>(flag); }

    inline constexpr void clear() { flag_ = 0; }
    
    friend std::ostream& operator<<(std::ostream& os, const StatusFlag& flag) { return os << flag.flag_; }

    int flag_;
};

namespace condition_detail{
    template <class T>
    constexpr inline T operator~(T a) { return (T) ~(int)a; }
    template <class T, class K>
    constexpr inline T operator|(T a, K b) { return (T)((int)a | (int)b); }
    template <class T, class K>
    constexpr inline T operator&(T a, K b) { return (T)((int)a & (int)b); }
    template <class T, class K>
    constexpr inline T operator^(T a, K b) { return (T)((int)a ^ (int)b); }
    template <class T, class K>
    constexpr inline T &operator|=(T &a, K b) { return (T &)((int &)a |= (int)b); }
    template <class T, class K>
    constexpr inline T &operator&=(T &a, K b) { return (T &)((int &)a &= (int)b); }
    template <class T, class K>
    constexpr inline T &operator^=(T &a, K b) { return (T &)((int &)a ^= (int)b); }

    template <class T, class K>
    bool not_set(const T flag, const K condition) { return int(flag & condition) == 0; }
    template <class T, class K>
    bool all_set(const T flag, const K condition) { return int(flag | ~condition) == ~0; }
    template <class T, class K>
    bool any_set(const T flag, const K condition) { return !condition || (int(flag) & int(condition)); }
 
}

struct ConditionFunctor
{
    ConditionFunctor() : not_conditions(0), and_conditions(0), or_conditions(0) {}                                                                                                                                                                    // Default, no conditions
    ConditionFunctor(int and_conditions, int not_conditions = 0, int or_conditions = 0) : not_conditions(not_conditions), and_conditions(and_conditions), or_conditions(or_conditions) {}                                                             // Custom conditions;
    ConditionFunctor(StatusFlag and_conditions, StatusFlag not_conditions = (StatusFlag)0, StatusFlag or_conditions = (StatusFlag)0) : not_conditions((int)not_conditions), and_conditions((int)and_conditions), or_conditions((int)or_conditions) {} // Custom conditions;
    ConditionFunctor(StatusEnum and_conditions, StatusEnum not_conditions = (StatusEnum)0, StatusEnum or_conditions = (StatusEnum)0) : not_conditions((int)not_conditions), and_conditions((int)and_conditions), or_conditions((int)or_conditions) {} // Custom conditions;
    
    template <size_t N, size_t M, size_t O>
    ConditionFunctor(std::array<int, N> and_conditions, std::array<int, M> not_conditions, std::array<int, O> or_conditions) : not_conditions(0), and_conditions(0), or_conditions(0)
    {
        for (int i : and_conditions)
            and_conditions |= i;
        for (int i : not_conditions)
            not_conditions |= i;
        for (int i : or_conditions)
            or_conditions |= i;
    }



    ConditionFunctor(const ConditionFunctor &other) = default;
    ConditionFunctor(ConditionFunctor &&other) = default;
    ConditionFunctor &operator=(const ConditionFunctor &other) = default;
    ConditionFunctor &operator=(ConditionFunctor &&other) = default;

    const int not_conditions;                                                                                                                                                                                                                          // Bitwise condtions that must not be met
    const int and_conditions;                                                                                                                                                                                                                          // Bitwise condtions that must be met
    const int or_conditions;                                                                                                                                                                                                                           // Bitwise condtions of which at least one must be met

    // First checks if the flags that must be met are met, then checks if the flags that must not be met are not met
    inline bool operator()(const StatusFlag flag) const
    {
        int iflag = (int)flag;
        return condition_detail::not_set(iflag, not_conditions) && condition_detail::all_set(iflag, and_conditions) && condition_detail::any_set(iflag, or_conditions);
    }
    inline constexpr bool operator()(const int flag) const { return condition_detail::not_set(flag, not_conditions) && condition_detail::all_set(flag, and_conditions) && condition_detail::any_set(flag, or_conditions); }
};