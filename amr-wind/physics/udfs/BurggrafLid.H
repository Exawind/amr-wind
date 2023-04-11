#ifndef BURGGRAF_LID_H
#define BURGGRAF_LID_H

#include "AMReX_Geometry.H"
#include "AMReX_Gpu.H"

namespace amr_wind {

class Field;

namespace udf {

struct BurggrafLid
{
    struct DeviceOp
    {
        int idir{0};

        AMREX_GPU_DEVICE
        inline void operator()(
            const amrex::IntVect& iv,
            amrex::Array4<amrex::Real> const& field,
            amrex::GeometryData const& geom,
            const amrex::Real /*unused*/,
            amrex::Orientation /*unused*/,
            const int comp,
            const int dcomp,
            const int orig_comp) const
        {
            const auto* problo = geom.ProbLo();
            const auto* dx = geom.CellSize();
            const auto x = problo[idir] + (iv[idir] + 0.5) * dx[idir];
            amrex::Real val =
                16 * (std::pow(x, 4) - 2 * std::pow(x, 3) + std::pow(x, 2));
            amrex::Real uvec[AMREX_SPACEDIM];
            uvec[0] = val;
            uvec[1] = 0.0;
            uvec[2] = 0.0;

            field(iv[0], iv[1], iv[2], dcomp + comp) = uvec[orig_comp + comp];
        }
    };
    using DeviceType = DeviceOp;

    static std::string identifier() { return "BurggrafLid"; }

    explicit BurggrafLid(const Field& fld);

    DeviceType device_instance() const { return m_op; }

    DeviceOp m_op;
};

} // namespace udf
} // namespace amr_wind

#endif /* BURGGRAF_LID_H */