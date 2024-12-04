#pragma once
#include <cassert>
#include <iterator>
#include <array>

template <typename T, size_t N>
struct MDSpan
{
    static_assert(N>0);
    using value_type = T;
    using pointer    = T*;
    using size_t  = std::size_t;
    using array_t = std::array<T,N>;

    inline constexpr MDSpan() : data_(nullptr) {}
    inline constexpr MDSpan(T *data, array_t shape) : data_(data), shape_(shape) {
            stride_[N-1] = 1;
            for(int axis=N-2;axis>=0;axis--) 
                stride_[axis] = stride_[axis+1]*shape_[axis+1];
    }
    inline constexpr MDSpan(T *begin, T *end):  MDSpan(begin, {end-begin}) {} 
    inline constexpr MDSpan(const MDSpan<T,N> &other) = default;
    inline constexpr MDSpan(MDSpan<T,N> &&other) = default;
    inline constexpr MDSpan(T& value) : MDSpan(&T,{1}) {}
    
    template <typename U>
    inline constexpr MDSpan<U,N> as_MDSpan() const {
        array_t shape = shape_;
        for(int axis=0;axis<N;axis++){ shape[i] *= sizeof(T); shape[i] /= sizeof(U); }
        return MDSpan<U,N>(reinterpret_cast<U*>(data_), shape);
    }
    inline constexpr size_t offset_of(array_t index) {
        size_t offset = 0;
#pragma unroll(3)
        for(int axis=0;axis<N;axis++) offset += index[axis]*stride_[axis];
        return offset;
    }
    inline constexpr array_t index_of(size_t offset) {
        array_t index;
        // Requirement: Either stride[i] = shape[0]   * ... * shape[i-1]
        //              or     stride[i] = shape[i+1] * ... * shape[n-1]
        if(stride_[0] < stride_[n-1]){
            for(int axis=N-1;axis>=0;axis--){
                index[i] = offset / stride[i];
                offset  %= stride[i];
            }
        } else { // stride_[0] > stride_[n-1]
            for(int axis=0;axis<N;axis++){
                index[i] = offset / stride[i];
                offset  %= stride[i];
            }
        }
    }
    
    inline constexpr MDSpan<T> subSpan(array_t index, array_t shape) const {
         size_t offset = offset_of(index);
         return MDSpan<T>(data_ + offset, shape); 
    }
    inline constexpr MDSpan<T,N>& operator= (const MDSpan<T,N> &other)  = default;
    // inline constexpr MDSpan<T,N>& operator= (const MDSpan<T,N> &other) { 
    //     data_ = other.data_; shape_ = other.shape_; stride_ = other.stride_; return *this; 
    // }
    inline constexpr MDSpan<T,N>& operator= (MDSpan<T,N> &&other) = default; // { return *this = other; }
    inline bool operator==(const MDSpan<T,N> other) {
        size_t size_ = size();
        for(size_t i=0;i<size_;i++){
            array_t index = index_of(i);
            other_.data_[index] == data_[index];
        }
    }

    inline constexpr T  operator[](const size_t index) const { return operator[]({ array_t{{index}} ); } 
    inline constexpr T& operator[](const size_t index)       { return operator[]({ array_t{{index}} ); }     

    inline constexpr T& operator[](const array_t &index) {
//#ifndef NDEBUG        
        for(int axis=0;axis<N;axis++) assert(index[i] < shape_[i]); // TODO: Langsomt, til debug naar virker
        assert(data_); 
//#endif
        return data_[offset_of(index)];
    }
    inline constexpr T operator[](const array_t &index) const {
//#ifndef NDEBUG        
        for(int axis=0;axis<N;axis++) assert(index[i] < shape_[i]); // TODO: Langsomt, til debug naar virker
        assert(data_); 
//#endif
        return data_[offset_of(index)];
    }

    inline constexpr T &at(const array_t &index) const {
        for(int axis=0;axis<N;axis++) assert(index[i] < shape_[i]); // Behold for evigt.
        assert(data_ != 0);         
        return data_[offset_of(index)];
    }
    inline constexpr T *data() const { return data_; }
    inline constexpr size_t size() const { 
        size_t size_ = 1;
        for(int axis=0;axis<N;axis++) size_ *= shape_[axis];
        return size_;
    }
    inline constexpr array_t shape() const { return shape_; }
    inline constexpr bool empty() const { return size() == 0; }
    inline constexpr T *begin() const { return data_; }
    inline constexpr T *end() const { return data_ + size(); }
    inline constexpr T &front() const { return data_[0]; }
    inline constexpr T &back() const { 
        array_t last_index = shape_;
        for(int axis=0;axis<N;axis++) last_index[axis]--;
        return data_[offset_of(last_index)]; 
    }
    template <typename U>
    friend std::ostream &operator<<(std::ostream &os, const MDSpan<U> &vec);

private:
    T *data_;
    array_t<size_t,N> shape_, strides_;
};
