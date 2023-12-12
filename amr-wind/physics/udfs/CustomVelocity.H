#ifndef CUSTOM_VELOCITY_H
#define CUSTOM_VELOCITY_H

#include "AMReX_Geometry.H"
#include "AMReX_Gpu.H"

namespace amr_wind {

class Field;

namespace udf {

struct CustomVelocity
{
    struct DeviceOp
    {

        AMREX_GPU_DEVICE
        inline void operator()(
            const amrex::IntVect& /*iv*/,
            amrex::Array4<amrex::Real> const& /*field*/,
            amrex::GeometryData const& /*geom*/,
            const amrex::Real /*time*/,
            amrex::Orientation /*ori*/,
            const int /*comp*/,
            const int /*dcomp*/,
            const int /*orig_comp*/) const
        {
            // Compute quantities to set the field values. For example:
            // clang-format off
            // const auto* problo = geom.ProbLo();
            // const auto* dx = geom.CellSize();
            // const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> vel = {1.0, 0.0, 0.0};
            // field(iv[0], iv[1], iv[2], dcomp + comp) = vel[orig_comp + comp];
            // clang-format on
        }
    };
    using DeviceType = DeviceOp;

    static std::string identifier() { return "CustomVelocity"; }

    explicit CustomVelocity(const Field& fld);

    DeviceType device_instance() const { return m_op; }

    DeviceOp m_op;
};

} // namespace udf
} // namespace amr_wind

#endif /* CUSTOM_VELOCITY_H */
