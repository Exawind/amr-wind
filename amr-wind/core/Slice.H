#ifndef SLICE_H
#define SLICE_H

#include <iterator>
#include "AMReX_Vector.H"
#include "AMReX_Gpu.H"

namespace amr_wind {
namespace utils {

template <typename T>
struct Slice
{
    T* ptr_begin;
    T* ptr_end;

    explicit Slice(T* pbegin, const size_t n)
        : ptr_begin(pbegin), ptr_end(pbegin + n)
    {}

    explicit Slice(T* pbegin, T* pend) : ptr_begin(pbegin), ptr_end(pend) {}

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T& operator[](const size_t idx)
    {
        return *(ptr_begin + idx);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const T&
    operator[](const size_t idx) const
    {
        return *(ptr_begin + idx);
    }

    size_t size() const { return std::distance(ptr_begin, ptr_end); }

    T* data() { return ptr_begin; }
    T* begin() { return ptr_begin; }
    T* end() { return ptr_end; }
    const T* data() const { return ptr_begin; }
    const T* begin() const { return ptr_begin; }
    const T* end() const { return ptr_end; }
    const T* cbegin() const { return ptr_begin; }
    const T* cend() const { return ptr_end; }
};

template <typename T>
inline Slice<T>
slice(std::vector<T>& vec, const size_t start, const size_t count)
{
    AMREX_ASSERT((start + count) <= vec.size());
    return Slice<T>{&vec[start], count};
}

template <typename T>
inline Slice<T> slice(std::vector<T>& vec, const size_t start)
{
    AMREX_ASSERT(start < vec.size());
    return slice(vec, start, vec.size() - start);
}

template <typename T>
inline Slice<T>
slice(amrex::Vector<T>& vec, const size_t start, const size_t count)
{
    AMREX_ASSERT((start + count) <= vec.size());
    return Slice<T>{&vec[start], count};
}

template <typename T>
inline Slice<T> slice(amrex::Vector<T>& vec, const size_t start)
{
    AMREX_ASSERT(start < vec.size());
    return slice(vec, start, vec.size() - start);
}

template <typename T>
inline Slice<T>
slice(amrex::Gpu::DeviceVector<T>& vec, const size_t start, const size_t count)
{
    AMREX_ASSERT((start + count) <= vec.size());
    return Slice<T>{(vec.data() + start), count};
}

template <typename T>
inline Slice<T> slice(amrex::Gpu::DeviceVector<T>& vec, const size_t start)
{
    AMREX_ASSERT(start < vec.size());
    return slice(vec, start, vec.size() - start);
}

} // namespace utils
} // namespace amr_wind

#endif /* SLICE_H */
