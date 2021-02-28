#ifndef LINEARPROFILE_H
#define LINEARPROFILE_H

#include "AMReX_Geometry.H"
#include "AMReX_Gpu.H"

namespace amr_wind {

class Field;

namespace udf {

struct LinearProfile
{
    struct DeviceOp
    {
        int idir{2};
        amrex::Real zmin;
        amrex::Real zmax;
        amrex::Real fld_min[AMREX_SPACEDIM];
        amrex::Real fld_max[AMREX_SPACEDIM];

        AMREX_GPU_DEVICE
        inline void operator()(
            const amrex::IntVect& iv,
            amrex::Array4<amrex::Real> const& field,
            amrex::GeometryData const& geom,
            const amrex::Real,
            amrex::Orientation,
            const int comp) const
        {
            const auto* problo = geom.ProbLo();
            const auto* dx = geom.CellSize();
            const auto zco = problo[idir] + (iv[idir] + 0.5) * dx[idir];

            amrex::Real val;
            if (zco < zmin) {
                val = fld_min[comp];
            } else if (zco > zmax) {
                val = fld_max[comp];
            } else {
                val = fld_min[comp] + (fld_max[comp] - fld_min[comp]) *
                                          (zco - zmin) / (zmax - zmin);
            }
            field(iv[0], iv[1], iv[2], comp) = val;
        }
    };

    using DeviceType = DeviceOp;

    static std::string identifier() { return "LinearProfile"; }

    explicit LinearProfile(const Field& fld);

    DeviceType device_instance() const { return m_op; }

    DeviceOp m_op;
};

} // namespace udf
} // namespace amr_wind

#endif /* LINEARPROFILE_H */
