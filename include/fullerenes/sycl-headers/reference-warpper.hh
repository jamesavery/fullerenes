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
    
    template <typename U>
    inline constexpr friend std::ostream& operator<<(std::ostream& os, const ReferenceWrapper<U>& ref);
};
