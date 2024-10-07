#pragma once
#include <cassert>

template <typename T>
struct SyclVector
{   
    using value_type = T;
    using pointer = T*;
    using size_t = std::size_t;
    //Constructors implementations depend on sycl, so they are not defined here
    SyclVector(size_t size);
    SyclVector(size_t size, T value);
    SyclVector(const SyclVector<T> &other);
    SyclVector<T> &operator=(const SyclVector<T> &other);
    ~SyclVector();

    void resize(size_t new_size);
    void resize(size_t new_size, T value);
    void resize(size_t new_size, size_t front, size_t back, size_t seg_size);

    void reserve(size_t new_capacity);

    //Movement semantics can be defined here
    SyclVector() : size_(0), capacity_(0), data_(nullptr) {}
    SyclVector(SyclVector<T> &&other) = default;
    SyclVector<T> &operator=(SyclVector<T> &&other) = default;

    inline constexpr operator Span<T>() const { return Span<T>(data_, size_); }
    inline constexpr void fill(T data) { std::fill(begin(), end(), data); }
    inline constexpr T *data() const { return data_; }
    inline constexpr size_t size() const { return size_; }
    inline constexpr size_t capacity() const { return capacity_; }
    
    
    inline constexpr void clear() { size_ = 0; }

    inline constexpr T &operator[](size_t index) { if(index >= size_) printf("Index: %d, Size: %d\n", index, size_); assert (index < size_); return data_[index]; }
    inline constexpr const T &operator[](size_t index) const { if (index >= size_) printf("Index: %d, Size: %d\n", index, size_); assert (index < size_); return data_[index]; }

    inline constexpr T &at(size_t index) { assert(index < size_); return data_[index]; }
    inline constexpr const T &at(size_t index) const { assert(index < size_); return data_[index]; }

    inline constexpr bool operator==(const SyclVector<T> &other) const {
        return Span<T>(*this) == Span<T>(other);
    }

    inline constexpr void push_back(const T &value) { 
        if(size_ == capacity_){
            size_t new_capacity = capacity_ == 0 ? 1 : 2*capacity_;
            reserve(new_capacity);
        }
        data_[size_++] = value;
    }


    inline constexpr T pop_back() { assert(size_ > 0); return data_[--size_]; }

    template <typename U>
    friend std::ostream &operator<<(std::ostream &os, const SyclVector<U> &vec);

    inline constexpr T *begin() const { return data_; }
    inline constexpr T *end() const { return data_ + size_; }

    inline constexpr T &back() const { return data_[size_ - 1]; }
    inline constexpr T &front() const { return data_[0]; }

private:
    size_t size_;
    size_t capacity_;
    pointer data_;
};