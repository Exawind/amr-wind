#ifndef VS_TENSOR_H
#define VS_TENSOR_H

#include "AMReX_Gpu.H"
#include "amr-wind/core/vs/vstraits.H"
#include "amr-wind/core/vs/vector.H"

namespace amr_wind::vs {

/** Tensor in 3D vector space
 */
template <typename T>
struct TensorT
{
    T vv[9]{Traits::zero(), Traits::zero(), Traits::zero(),
            Traits::zero(), Traits::zero(), Traits::zero(),
            Traits::zero(), Traits::zero(), Traits::zero()};

    static constexpr int ncomp = 9;
    using size_type = int;
    using value_type = T;
    using reference = T&;
    using iterator = T*;
    using const_iterator = const T*;
    using Traits = DTraits<T>;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr TensorT() = default;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr TensorT(
        const T& xx,
        const T& xy,
        const T& xz,
        const T& yx,
        const T& yy,
        const T& yz,
        const T& zx,
        const T& zy,
        const T& zz)
        : vv{xx, xy, xz, yx, yy, yz, zx, zy, zz}
    {}

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE TensorT(
        const VectorT<T>& x,
        const VectorT<T>& y,
        const VectorT<T>& z,
        const bool transpose = false);

    ~TensorT() = default;
    TensorT(const TensorT&) = default;
    TensorT(TensorT&&) = default;
    TensorT& operator=(const TensorT&) & = default;
    TensorT& operator=(const TensorT&) && = delete;
    TensorT& operator=(TensorT&&) & = default;
    TensorT& operator=(TensorT&&) && = delete;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static constexpr TensorT<T>
    zero() noexcept
    {
        return TensorT<T>{Traits::zero(), Traits::zero(), Traits::zero(),
                          Traits::zero(), Traits::zero(), Traits::zero(),
                          Traits::zero(), Traits::zero(), Traits::zero()};
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static constexpr TensorT<T>
    I() noexcept
    {
        return TensorT{Traits::one(),  Traits::zero(), Traits::zero(),
                       Traits::zero(), Traits::one(),  Traits::zero(),
                       Traits::zero(), Traits::zero(), Traits::one()};
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void rows(
        const VectorT<T>& x, const VectorT<T>& y, const VectorT<T>& z) noexcept;
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void cols(
        const VectorT<T>& x, const VectorT<T>& y, const VectorT<T>& z) noexcept;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T> x() const noexcept;
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T> y() const noexcept;
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T> z() const noexcept;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T> cx() const noexcept;
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T> cy() const noexcept;
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T> cz() const noexcept;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T& xx() & noexcept
    {
        return vv[0];
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T& xy() & noexcept
    {
        return vv[1];
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T& xz() & noexcept
    {
        return vv[2];
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T& yx() & noexcept
    {
        return vv[3];
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T& yy() & noexcept
    {
        return vv[4];
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T& yz() & noexcept
    {
        return vv[5];
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T& zx() & noexcept
    {
        return vv[6];
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T& zy() & noexcept
    {
        return vv[7];
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T& zz() & noexcept
    {
        return vv[8];
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const T& xx() const& noexcept
    {
        return vv[0];
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const T& xy() const& noexcept
    {
        return vv[1];
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const T& xz() const& noexcept
    {
        return vv[2];
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const T& yx() const& noexcept
    {
        return vv[3];
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const T& yy() const& noexcept
    {
        return vv[4];
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const T& yz() const& noexcept
    {
        return vv[5];
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const T& zx() const& noexcept
    {
        return vv[6];
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const T& zy() const& noexcept
    {
        return vv[7];
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const T& zz() const& noexcept
    {
        return vv[8];
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T& operator[](size_type pos) &
    {
        return vv[pos];
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const T&
    operator[](size_type pos) const&
    {
        return vv[pos];
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T* data() noexcept { return vv; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const T* data() const noexcept
    {
        return vv;
    }

    iterator begin() noexcept { return vv; }
    iterator end() noexcept { return vv + ncomp; }
    const_iterator cbegin() const noexcept { return vv; }
    const_iterator cend() const noexcept { return vv + ncomp; }
    size_type size() const noexcept { return ncomp; }
};

using Tensor = TensorT<amrex::Real>;

} // namespace amr_wind::vs

#include "amr-wind/core/vs/tensorI.H"

#endif /* VS_TENSOR_H */
