#include "amr-wind/physics/udfs/Rankine.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/core/vs/vector.H"
#include "amr-wind/equation_systems/icns/icns.H"

#include "AMReX_ParmParse.H"

namespace amr_wind::udf {

Rankine::Rankine(const Field& fld)
{
    AMREX_ALWAYS_ASSERT(fld.name() == pde::ICNS::var_name());
    AMREX_ALWAYS_ASSERT(fld.num_comp() == AMREX_SPACEDIM);

    const int ncomp = fld.num_comp();

    {
        amrex::ParmParse pp("incflo");
        amrex::Vector<amrex::Real> vel(0.0, ncomp);
        pp.getarr("velocity", vel);
        AMREX_ALWAYS_ASSERT(vel.size() == ncomp);
        for (int i = 0; i < ncomp; ++i) {
            m_op.vel_ref[i] = vel[i];
        }
    }
    {
        amrex::ParmParse pp("Rankine");
        pp.query("Umax", m_op.Umax);
        pp.query("Rmax", m_op.Rmax);
    }
}

} // namespace amr_wind::udf
