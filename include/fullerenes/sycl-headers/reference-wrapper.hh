#pragma once
#include <functional>
#include <iostream>

template<typename T>
struct ReferenceWrapper : public std::reference_wrapper<T> {
    using std::reference_wrapper<T>::reference_wrapper;
    using std::reference_wrapper<T>::get;

    ReferenceWrapper() = delete;
    ReferenceWrapper(const ReferenceWrapper<T>& other) = default;
    ReferenceWrapper(const std::reference_wrapper<T>& other) : std::reference_wrapper<T>(other) {}
    ReferenceWrapper(const T& other) : std::reference_wrapper<T>(other) {}
    ReferenceWrapper<T>& operator=(const ReferenceWrapper<T>& other) = default;
    ReferenceWrapper<T>& operator=(const T& other) { get() = other; return *this; }



    ReferenceWrapper<T>& operator=(const std::reference_wrapper<T>& other) { get() = other.get(); return *this; }
    ReferenceWrapper<T>& operator=(std::reference_wrapper<T>&& other) { get() = std::move(other.get()); return *this; }

    T* operator->() const noexcept { return &(get()); }
    T& operator*() const noexcept { return get(); }
    T& operator[](int index) const { return get()[index]; }
    


    bool operator ==(const ReferenceWrapper<T>& other) const { return get() == other.get(); }
    bool operator !=(const ReferenceWrapper<T>& other) const { return get() != other.get(); }
    bool operator <(const ReferenceWrapper<T>& other) const { return get() < other.get(); }
    bool operator >(const ReferenceWrapper<T>& other) const { return get() > other.get(); }
    bool operator <=(const ReferenceWrapper<T>& other) const { return get() <= other.get(); }
    bool operator >=(const ReferenceWrapper<T>& other) const { return get() >= other.get(); }

    bool operator ==(const T& other) const { return get() == other; }
    bool operator !=(const T& other) const { return get() != other; }
    bool operator <(const T& other) const { return get() < other; }
    bool operator >(const T& other) const { return get() > other; }
    bool operator <=(const T& other) const { return get() <= other; }
    bool operator >=(const T& other) const { return get() >= other; }

    bool operator ==(const std::reference_wrapper<T>& other) const { return get() == other.get(); }
    bool operator !=(const std::reference_wrapper<T>& other) const { return get() != other.get(); }
    bool operator <(const std::reference_wrapper<T>& other) const { return get() < other.get(); }
    bool operator >(const std::reference_wrapper<T>& other) const { return get() > other.get(); }
    bool operator <=(const std::reference_wrapper<T>& other) const { return get() <= other.get(); }
    bool operator >=(const std::reference_wrapper<T>& other) const { return get() >= other.get(); }
    
    // Host-only, like its definition: std::ostream is not a literal type, so
    // this can never be constant-evaluated, and under clang CUDA a constexpr
    // function is implicitly __host__ __device__ -- which would both drag
    // std::ostream onto the device and make this friend declaration disagree
    // with the definition.
    template <typename U>
    inline friend std::ostream& operator<<(std::ostream& os, const ReferenceWrapper<U>& ref);
};
