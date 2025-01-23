#pragma once
#include <cassert>
#include <iterator>
#include <array>

template <typename T>
struct is_std_array : std::false_type {};

template <typename T, std::size_t N>
struct is_std_array<std::array<T, N>> : std::true_type {};

template <typename T>
struct Span
{   
    using value_type = T;
    using pointer = T*;
    using size_t = std::size_t;
    inline constexpr Span() : data_(nullptr), size_(0) {}
    inline constexpr Span(T *data, size_t size) : data_(data), size_(size) {}
    inline constexpr Span(T *begin, T *end) : data_(begin), size_(std::distance(begin, end)) {}
    inline constexpr Span(const Span<T> &other) = default;
    inline constexpr Span(Span<T> &&other) = default;
    inline constexpr Span(const Span<T>& other, size_t offset) : data_(other.data_ + offset), size_(other.size_ - offset) {}
    inline constexpr Span(Span<T>&& other, size_t offset) : data_(other.data_ + offset), size_(other.size_ - offset) {}
    inline constexpr Span(T& value) : data_(&value), size_(1) {}
    
    template <typename U>
    inline constexpr Span<U> as_span() const {
        return Span<U>(reinterpret_cast<U*>(data_), (sizeof(T) * size_ / sizeof(U)) );
    }

    inline constexpr Span<T> subspan(size_t offset) const { assert(size_ - offset >= 0); return Span<T>(data_ + offset, size_ - offset); }
    inline constexpr Span<T> subspan(size_t offset, size_t count) const { assert(offset + count <= size_);  return Span<T>(data_ + offset, count); }
    inline constexpr Span<T>& operator= (const Span<T> &other) { data_ = other.data_; size_ = other.size_; return *this; }
    inline constexpr Span<T>& operator= (Span<T> &&other) { return *this = other; }
    inline bool operator==(const Span<T> other) const;
    inline constexpr T &operator[](size_t index) const {assert(index < size_); assert(data_); return data_[index];}
    inline constexpr T &at(size_t index) const{assert(index < size_); assert(data_); return data_[index];}
    inline constexpr T *data() const { return data_; }
    inline constexpr size_t size() const { return size_; }
    inline constexpr bool empty() const { return size_ == 0; }
    inline constexpr size_t size_bytes() const { return size_ * sizeof(T); }
    inline constexpr T *begin() const { return data_; }
    inline constexpr T *end() const { return data_ + size_; }
    inline constexpr T &front() const { return data_[0]; }
    inline constexpr T &back() const { return data_[size_ - 1]; }
    template <typename U>
    friend std::ostream &operator<<(std::ostream &os, const Span<U> &vec);
    
private:
    T *data_;
    size_t size_;
};

template <typename T>
Span(T*, typename Span<T>::size_t) -> Span<T>;

template <typename T>
Span(T*, std::size_t) -> Span<T>;

template <typename T>
Span(T*, T*) -> Span<T>;

template <typename T>
Span(T&) -> Span<T>;

