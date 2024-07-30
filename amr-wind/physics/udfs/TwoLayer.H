#ifndef TWO_LAYER_H
#define TWO_LAYER_H

#include "AMReX_Geometry.H"
#include "AMReX_Gpu.H"

namespace amr_wind {

class Field;

namespace udf {

struct TwoLayer
{
    struct DeviceOp
    {
        // velocities of the top and bottom layer respectively
        //amrex::Real pvel{0.0};  amrex::Real mvel{0.0};
        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> pvel = {0.0};
        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> mvel = {0.0};

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
            const auto* problo = geom.ProbLo();
            const auto* dx = geom.CellSize();
            const auto z = problo[2] + (iv[2] + 0.5) * dx[2];   // not true for w-mac vels

            amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> vel;
            if (z>=0.5) {
                vel = {pvel[0], pvel[1], 0.0};
            } else {
                vel = {mvel[0], mvel[1], 0.0};
            }

            // perturb by some factor at time 0
            if (time == 0.0)  {
                vel[0] *= 0.9;
                vel[1] *= 0.9;
            }

            field(iv[0], iv[1], iv[2], dcomp + comp) = vel[orig_comp + comp];
        }
    };
    using DeviceType = DeviceOp;

    static std::string identifier() { return "TwoLayer"; }

    explicit TwoLayer(const Field& fld);

    DeviceType device_instance() const { return m_op; }

    DeviceOp m_op;
};

} // namespace udf
} // namespace amr_wind

#endif /* TWO_LAYER_H */
