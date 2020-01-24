#include <incflo.H>

using namespace amrex;

DiffusionTensorOp*
incflo::get_diffusion_tensor_op ()
{
    if (!m_diffusion_tensor_op) m_diffusion_tensor_op.reset(new DiffusionTensorOp(this));
    return m_diffusion_tensor_op.get();
}

DiffusionScalarOp*
incflo::get_diffusion_scalar_op ()
{
    if (!m_diffusion_scalar_op) m_diffusion_scalar_op.reset(new DiffusionScalarOp(this));
    return m_diffusion_scalar_op.get();
}

Vector<Array<LinOpBCType,AMREX_SPACEDIM> >
incflo::get_diffuse_tensor_bc (Orientation::Side side) const noexcept
{
    Array<LinOpBCType,AMREX_SPACEDIM> r;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (Geom(0).isPeriodic(dir)) {
            r[dir] = LinOpBCType::Periodic;
        } else {
            auto bc = m_bc_type[Orientation(dir,side)];
            switch (bc)
            {
            case BC::pressure_inflow:
            case BC::pressure_outflow:
            {
                r[dir] = LinOpBCType::Neumann;
                break;
            }
            case BC::mass_inflow:
            case BC::no_slip_wall:
            {
                r[dir] = LinOpBCType::Dirichlet;
                break;
            }
            default:
                amrex::Abort("get_diffuse_tensor_bc: undefined BC type");
            };
        }
    }
    return {r, r, r};
}

Array<LinOpBCType,AMREX_SPACEDIM>
incflo::get_diffuse_scalar_bc (Orientation::Side side) const noexcept
{
    Array<LinOpBCType,AMREX_SPACEDIM> r;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (Geom(0).isPeriodic(dir)) {
            r[dir] = LinOpBCType::Periodic;
        } else {
            auto bc = m_bc_type[Orientation(dir,side)];
            switch (bc)
            {
            case BC::pressure_inflow:
            case BC::pressure_outflow:
            case BC::no_slip_wall:
            {
                r[dir] = LinOpBCType::Neumann;
                break;
            }
            case BC::mass_inflow:
            {
                r[dir] = LinOpBCType::Dirichlet;
                break;
            }
            default:
                amrex::Abort("get_diffuse_tensor_bc: undefined BC type");
            };
        }
    }
    return r;
}

Array<MultiFab,AMREX_SPACEDIM>
incflo::average_velocity_eta_to_faces (int lev, MultiFab const& cc_eta)
{
    const auto& ba = cc_eta.boxArray();
    const auto& dm = cc_eta.DistributionMap();
    const auto& fact = cc_eta.Factory();
    Array<MultiFab,AMREX_SPACEDIM> r{MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(0)),
                                              dm, 1, 0, MFInfo(), fact),
                                     MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(1)),
                                              dm, 1, 0, MFInfo(), fact),
                                     MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(2)),
                                              dm, 1, 0, MFInfo(), fact)};
    amrex::average_cellcenter_to_face(GetArrOfPtrs(r), cc_eta, Geom(lev));
    return r;
}

Array<MultiFab,AMREX_SPACEDIM>
incflo::average_tracer_eta_to_faces (int lev, int comp, MultiFab const& cc_eta)
{
    const auto& ba = cc_eta.boxArray();
    const auto& dm = cc_eta.DistributionMap();
    const auto& fact = cc_eta.Factory();
    MultiFab cc(cc_eta, amrex::make_alias, comp, 1);
    Array<MultiFab,AMREX_SPACEDIM> r{MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(0)),
                                              dm, 1, 0, MFInfo(), fact),
                                     MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(1)),
                                              dm, 1, 0, MFInfo(), fact),
                                     MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(2)),
                                              dm, 1, 0, MFInfo(), fact)};
    amrex::average_cellcenter_to_face(GetArrOfPtrs(r), cc_eta, Geom(lev));
    return r;
}

