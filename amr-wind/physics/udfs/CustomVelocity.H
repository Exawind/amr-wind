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
        // clang-format off
        // Declare parameters here if needed. For example:
        // amrex::Real foo{1.0};
        amrex::Real pvel{1.0};  amrex::Real mvel{0.5};
        // amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> bar = {0.0};
        // clang-format on

        AMREX_GPU_DEVICE
        inline void operator()(
            const amrex::IntVect& iv,
            amrex::Array4<amrex::Real> const& field,
            amrex::GeometryData const& geom,
            const amrex::Real time,
            amrex::Orientation ori,
            const int comp,
            const int dcomp,
            const int orig_comp) const
        {
            if ((iv[0] == -1) && (iv[1] == -1) && (iv[2] == 0))
                amrex::Print() << "****** filling custom velocity at " << iv
                               << " and at other ghost cells... " << std::endl;
            // Compute quantities to set the field values. For example:
            // clang-format off
            const auto* problo = geom.ProbLo();
            const auto* dx = geom.CellSize();
            const auto x = problo[0] + (iv[0] + 0.5) * dx[0];
            // const auto y = problo[1] + (iv[1] + 0.5) * dx[1];
            const auto z = problo[2] + (iv[2] + 0.5) * dx[2];
            // const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> vel = {1.0, 0.0, 0.0};

            // Once the above is done, fill the field as:
            amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> vel;
            if (x <= 0.5) {
                // at the low-x boundary, +ve vel in top portion, garbage in the rest
                //vel = {((z >= 0.5) ? pvel : -1e10), 0.0, 0.0};
                //vel = {((z >= 0.5) ? pvel : 0.0), 0.0, 0.0};
                vel = {((z >= 0.5) ? pvel : mvel), 0.0, 0.0};
            } else {
                // at the high-x boundary, -ve vel in bottom portion, garbage in the rest
                //vel = {((z < 0.5) ? mvel : 1e10), 0.0, 0.0};
                //vel = {((z < 0.5) ? mvel : 0.0), 0.0, 0.0};
                vel = {((z < 0.5) ? mvel : pvel), 0.0, 0.0};
            }

            field(iv[0], iv[1], iv[2], dcomp + comp) = vel[orig_comp + comp];
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
