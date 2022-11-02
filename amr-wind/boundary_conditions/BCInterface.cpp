#include "amr-wind/boundary_conditions/BCInterface.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/boundary_conditions/FixedGradientBC.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

BCIface::BCIface(Field& field) : m_field(field) {}

void BCIface::operator()(const amrex::Real value)
{
    amrex::Print() << "Initializing boundary conditions for " << m_field.name()
                   << std::endl;
    set_default_value(value);
    read_bctype();
    set_bcrec();
    read_values();
    set_bcfuncs();
    m_field.copy_bc_to_device();
}

inline void BCIface::set_default_value(const amrex::Real value)
{
    auto& bcval = m_field.bc_values();
    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        for (int i = 0; i < m_field.num_comp(); ++i) {
            bcval[ori][i] = value;
        }
    }
}

void BCIface::read_bctype()
{
    const std::string key = m_field.name() + "_type";
    auto& ibctype = m_field.bc_type();
    const auto& geom = m_field.repo().mesh().Geom(0);
    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        const auto& bcid = bcnames[ori];
        amrex::ParmParse pp(bcid);
        std::string bcstr = "null";
        pp.query("type", bcstr);
        pp.query(key.c_str(), bcstr);
        bcstr = amrex::toLower(bcstr);

        // Protect against copy/paste errors where user intended to add a BC but
        // forgot to turn off periodic in that direction, or vice versa.
        if (geom.isPeriodic(ori.coordDir()) && bcstr != "null") {
            amrex::Abort(
                "Setting " + bcstr + " BC on a periodic face " + bcid +
                " is not allowed");
        }

        // Process and quit early if this face is periodic
        if (geom.isPeriodic(ori.coordDir())) {
            ibctype[ori] = BC::periodic;
            continue;
        }

        if ((bcstr == "pressure_inflow") || (bcstr == "pi")) {
            ibctype[ori] = BC::pressure_inflow;
        } else if ((bcstr == "pressure_outflow") || (bcstr == "po")) {
            ibctype[ori] = BC::pressure_outflow;
        } else if ((bcstr == "mass_inflow") || (bcstr == "mi")) {
            ibctype[ori] = BC::mass_inflow;
        } else if ((bcstr == "no_slip_wall") || (bcstr == "nsw")) {
            ibctype[ori] = BC::no_slip_wall;
        } else if ((bcstr == "slip_wall") || (bcstr == "sw")) {
            ibctype[ori] = BC::slip_wall;
        } else if ((bcstr == "wall_model") || (bcstr == "wm")) {
            ibctype[ori] = BC::wall_model;
        } else if ((bcstr == "zero_gradient") || (bcstr == "zg")) {
            ibctype[ori] = BC::zero_gradient;
        } else if ((bcstr == "fixed_gradient") || (bcstr == "fg")) {
            ibctype[ori] = BC::fixed_gradient;
        } else {
            ibctype[ori] = BC::undefined;
        }

        if (ibctype[ori] == BC::undefined) {
            amrex::Abort(
                "Invalid BC specification for non-periodic boundary = " +
                bcstr);
        }
    }
}

void BCIface::set_bcfuncs()
{
    const auto& ibctype = m_field.bc_type();
    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        const auto bct = ibctype[ori];

        switch (bct) {
        case BC::fixed_gradient:
            m_field.register_custom_bc<FixedGradientBC>(ori);
            break;

        default:
            break;
        }
    }
}

std::pair<const std::string, const std::string> BCIface::get_dirichlet_udfs()
{
    const auto& fname = m_field.name();
    const auto& bctype = m_field.bc_type();
    const std::string inflow_key = fname + ".inflow_type";
    const std::string wall_key = fname + ".wall_type";
    std::string inflow_udf{"ConstDirichlet"};
    std::string wall_udf{"ConstDirichlet"};
    bool has_inflow_udf = false;
    bool has_wall_udf = false;

    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        const auto& bcid = bcnames[ori];
        const auto bct = bctype[ori];
        amrex::ParmParse pp(bcid);

        switch (bct) {
        case BC::mass_inflow: {
            if (pp.contains(inflow_key.c_str())) {
                std::string val;
                pp.get(inflow_key.c_str(), val);

                if (has_inflow_udf && (inflow_udf != val)) {
                    amrex::Abort(
                        "BCVelocity: Inflow UDF must be same for all inflow "
                        "faces");
                } else {
                    inflow_udf = val;
                }
            }
            break;
        }

        case BC::slip_wall: {
            if (pp.contains(wall_key.c_str())) {
                std::string val;
                pp.get(wall_key.c_str(), val);

                if (has_wall_udf && (wall_udf != val)) {
                    amrex::Abort(
                        "BCVelocity: Wall UDF must be same for all wall "
                        "faces");
                } else {
                    wall_udf = val;
                }
            }
            break;
        }

        default:
            break;
        }
    }

    return {inflow_udf, wall_udf};
}

void BCVelocity::set_bcrec()
{
    const auto& ibctype = m_field.bc_type();
    auto& bcrec = m_field.bcrec();

    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        const auto side = ori.faceDir();
        const auto bct = ibctype[ori];
        const int dir = ori.coordDir();

        switch (bct) {
        case BC::periodic:
            if (side == amrex::Orientation::low) {
                set_bcrec_lo(dir, amrex::BCType::int_dir);
            } else {
                set_bcrec_hi(dir, amrex::BCType::int_dir);
            }
            break;

        case BC::pressure_inflow:
        case BC::pressure_outflow:
        case BC::zero_gradient:
            if (side == amrex::Orientation::low) {
                set_bcrec_lo(dir, amrex::BCType::foextrap);
            } else {
                set_bcrec_hi(dir, amrex::BCType::foextrap);
            }
            break;

        case BC::mass_inflow:
        case BC::no_slip_wall:
            if (side == amrex::Orientation::low) {
                set_bcrec_lo(dir, amrex::BCType::ext_dir);
            } else {
                set_bcrec_hi(dir, amrex::BCType::ext_dir);
            }
            break;

        case BC::slip_wall:
        case BC::wall_model:
            if (side == amrex::Orientation::low) {
                // Tangential directions use higher-order extrapolation
                set_bcrec_lo(dir, amrex::BCType::hoextrap);
                // Normal direction uses dirichlet (ext-dir)
                bcrec[dir].setLo(dir, amrex::BCType::ext_dir);
            } else {
                // Tangential directions use higher-order extrapolation
                set_bcrec_hi(dir, amrex::BCType::hoextrap);
                // Normal direction uses dirichlet (ext-dir)
                bcrec[dir].setHi(dir, amrex::BCType::ext_dir);
            }
            break;

        default:
            amrex::Abort("Invalid incflo BC type encountered");
        }
    }
}

void BCVelocity::read_values()
{
    const auto& fname = m_field.name();
    const auto& bctype = m_field.bc_type();
    auto& bcval = m_field.bc_values();
    const int ndim = m_field.num_comp();
    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        const auto& bcid = bcnames[ori];
        const auto bct = bctype[ori];

        amrex::ParmParse pp(bcid);
        switch (bct) {
        case BC::mass_inflow:
            pp.queryarr(fname.c_str(), bcval[ori], 0, ndim);
            break;

        case BC::no_slip_wall:
            pp.queryarr(fname.c_str(), bcval[ori], 0, ndim);
            // Set normal component to zero
            bcval[ori][ori.coordDir()] = 0.0;
            break;

        default:
            pp.queryarr(fname.c_str(), bcval[ori], 0, ndim);
            break;
        }
    }
}

void BCScalar::set_bcrec()
{
    const auto& ibctype = m_field.bc_type();
    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        const auto side = ori.faceDir();
        const auto bct = ibctype[ori];
        const int dir = ori.coordDir();

        switch (bct) {
        case BC::periodic:
            if (side == amrex::Orientation::low) {
                set_bcrec_lo(dir, amrex::BCType::int_dir);
            } else {
                set_bcrec_hi(dir, amrex::BCType::int_dir);
            }
            break;

        case BC::pressure_inflow:
        case BC::pressure_outflow:
        case BC::zero_gradient:
            if (side == amrex::Orientation::low) {
                set_bcrec_lo(dir, amrex::BCType::foextrap);
            } else {
                set_bcrec_hi(dir, amrex::BCType::foextrap);
            }
            break;

        case BC::mass_inflow:
        case BC::no_slip_wall:
            if (side == amrex::Orientation::low) {
                set_bcrec_lo(dir, amrex::BCType::ext_dir);
            } else {
                set_bcrec_hi(dir, amrex::BCType::ext_dir);
            }
            break;

        case BC::slip_wall:
        case BC::wall_model:
        case BC::fixed_gradient:
            if (side == amrex::Orientation::low) {
                set_bcrec_lo(dir, amrex::BCType::hoextrap);
            } else {
                set_bcrec_hi(dir, amrex::BCType::hoextrap);
            }
            break;

        default:
            amrex::Abort("Invalid incflo BC type encountered");
        }
    }
}

void BCScalar::read_values()
{
    const auto& fname = m_field.name();
    const auto& bctype = m_field.bc_type();
    auto& bcval = m_field.bc_values();
    const int ndim = m_field.num_comp();
    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        const auto& bcid = bcnames[ori];
        const auto bct = bctype[ori];

        amrex::ParmParse pp(bcid);
        switch (bct) {
        case BC::mass_inflow:
            pp.getarr(fname.c_str(), bcval[ori], 0, ndim);
            break;

        default:
            pp.queryarr(fname.c_str(), bcval[ori], 0, ndim);
            break;
        }
    }
}

void BCPressure::read_values()
{
    const auto& fname = m_field.name();
    const auto& bctype = m_field.bc_type();
    auto& bcval = m_field.bc_values();
    const int ndim = m_field.num_comp();
    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        const auto& bcid = bcnames[ori];
        const auto bct = bctype[ori];

        amrex::ParmParse pp(bcid);
        switch (bct) {
        case BC::pressure_inflow:
        case BC::pressure_outflow:
            pp.queryarr(fname.c_str(), bcval[ori], 0, ndim);
            break;

        default:
            break;
        }
    }
}

void BCSrcTerm::set_bcrec()
{
    const auto& ibctype = m_field.bc_type();
    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        const auto side = ori.faceDir();
        const auto bct = ibctype[ori];
        const int dir = ori.coordDir();

        switch (bct) {
        case BC::periodic:
            if (side == amrex::Orientation::low) {
                set_bcrec_lo(dir, amrex::BCType::int_dir);
            } else {
                set_bcrec_hi(dir, amrex::BCType::int_dir);
            }
            break;

        default:
            if (side == amrex::Orientation::low) {
                set_bcrec_lo(dir, amrex::BCType::foextrap);
            } else {
                set_bcrec_hi(dir, amrex::BCType::foextrap);
            }
            break;
        }
    }
}

void BCFillPatchExtrap::set_bcrec()
{
    const auto& ibctype = m_field.bc_type();
    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        const auto side = ori.faceDir();
        const auto bct = ibctype[ori];
        const int dir = ori.coordDir();

        switch (bct) {
        case BC::periodic:
            if (side == amrex::Orientation::low) {
                set_bcrec_lo(dir, amrex::BCType::int_dir);
            } else {
                set_bcrec_hi(dir, amrex::BCType::int_dir);
            }
            break;

        default:
            if (side == amrex::Orientation::low) {
                set_bcrec_lo(dir, m_extrap_type);
            } else {
                set_bcrec_hi(dir, m_extrap_type);
            }
            break;
        }
    }
}

} // namespace amr_wind
