#include <incflo.H>

using namespace amrex;

DiffusionTensorOp*
incflo::get_diffusion_tensor_op ()
{
    if (!diffusion_tensor_op) diffusion_tensor_op.reset(new DiffusionTensorOp(this));
    return diffusion_tensor_op.get();
}

Vector<Array<LinOpBCType,AMREX_SPACEDIM> >
incflo::get_diffuse_tensor_bc (Orientation::Side side) const noexcept
{
    Vector<Array<LinOpBCType,AMREX_SPACEDIM> > r(AMREX_SPACEDIM);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (Geom(0).isPeriodic(dir)) {
            for (int vcomp = 0; vcomp < AMREX_SPACEDIM; ++vcomp) {
                r[vcomp][dir] = LinOpBCType::Periodic;
            }
        } else {
            auto bc = m_bc_type[Orientation(dir,side)];
            switch (bc)
            {
            case BC::pressure_inflow:
            case BC::pressure_outflow:
            {
                for (int vcomp = 0; vcomp < AMREX_SPACEDIM; ++vcomp) {
                    r[vcomp][dir] = LinOpBCType::Neumann;
                }        
                break;
            }
            case BC::mass_inflow:
            case BC::no_slip_wall:
            {
                for (int vcomp = 0; vcomp < AMREX_SPACEDIM; ++vcomp) {
                    r[vcomp][dir] = LinOpBCType::Dirichlet;
                }        
                break;
            }
            default:
                amrex::Abort("get_diffuse_tensor_bc: undefined BC type");
            };
        }
    }
    return r;
}
