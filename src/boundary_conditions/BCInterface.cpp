#include "BCInterface.H"
#include "FieldRepo.H"
#include "FixedGradientBC.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace {
amrex::Vector<std::string> bcnames = {"xlo", "ylo", "zlo", "xhi", "yhi", "zhi"};
}

BCIface::BCIface(Field& field)
    : m_field(field)
{}

inline void BCIface::set_bcrec_lo(
    int dir, amrex::BCType::mathematicalBndryTypes bcrec)
{
    auto& fbcrec = m_field.bcrec();
    for (int i = 0; i < m_field.num_comp(); ++i) fbcrec[i].setLo(dir, bcrec);
}

inline void BCIface::set_bcrec_hi(
    int dir, amrex::BCType::mathematicalBndryTypes bcrec)
{
    auto& fbcrec = m_field.bcrec();
    for (int i = 0; i < m_field.num_comp(); ++i) fbcrec[i].setHi(dir, bcrec);
}

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
    for (amrex::OrientationIter oit; oit; ++oit) {
        auto ori = oit();
        for (int i = 0; i < m_field.num_comp(); ++i) bcval[ori][i] = value;
    }
}

void BCIface::read_bctype()
{
    const std::string key = m_field.name() + "_type";
    auto& ibctype = m_field.bc_type();
    auto& geom = m_field.repo().mesh().Geom(0);
    for (amrex::OrientationIter oit; oit; ++oit) {
        auto ori = oit();
        // Process and quit early if this face is periodic
        if (geom.isPeriodic(ori.coordDir())) {
            ibctype[ori] = BC::periodic;
            continue;
        }

        const auto& bcid = bcnames[ori];
        amrex::ParmParse pp(bcid);
        std::string bcstr = "null";
        pp.query("type", bcstr);
        pp.query(key.c_str(), bcstr);
        bcstr = amrex::toLower(bcstr);

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

        if (ibctype[ori] == BC::undefined)  {
            amrex::Abort("No BC specified for non-periodic boundary");
        }
    }
}

void BCIface::set_bcfuncs()
{
    auto& ibctype = m_field.bc_type();
    for (amrex::OrientationIter oit; oit; ++oit) {
        auto ori = oit();
        const auto bct = ibctype[ori];

        switch(bct) {
        case BC::fixed_gradient:
            m_field.register_custom_bc<FixedGradientBC>(ori);
            break;

        default:
            break;
        }
    }
}

void BCVelocity::set_bcrec()
{
    auto& ibctype = m_field.bc_type();
    auto& bcrec = m_field.bcrec();

    for (amrex::OrientationIter oit; oit; ++oit) {
        auto ori = oit();
        const auto side = ori.faceDir();
        const auto bct = ibctype[ori];
        const int dir = ori.coordDir();

        switch (bct) {
        case BC::periodic:
            if (side == amrex::Orientation::low)
                set_bcrec_lo(dir, amrex::BCType::int_dir);
            else
                set_bcrec_hi(dir, amrex::BCType::int_dir);
            break;

        case BC::pressure_inflow:
        case BC::pressure_outflow:
        case BC::zero_gradient:
            if (side == amrex::Orientation::low)
                set_bcrec_lo(dir, amrex::BCType::foextrap);
            else
                set_bcrec_hi(dir, amrex::BCType::foextrap);
            break;

        case BC::mass_inflow:
        case BC::no_slip_wall:
            if (side == amrex::Orientation::low)
                set_bcrec_lo(dir, amrex::BCType::ext_dir);
            else
                set_bcrec_hi(dir, amrex::BCType::ext_dir);
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
    auto& fname = m_field.name();
    auto& bctype = m_field.bc_type();
    auto& bcval = m_field.bc_values();
    const int ndim = m_field.num_comp();
    for (amrex::OrientationIter oit; oit; ++oit) {
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
    auto& ibctype = m_field.bc_type();
    for (amrex::OrientationIter oit; oit; ++oit) {
        auto ori = oit();
        const auto side = ori.faceDir();
        const auto bct = ibctype[ori];
        const int dir = ori.coordDir();

        switch (bct) {
        case BC::periodic:
            if (side == amrex::Orientation::low)
                set_bcrec_lo(dir, amrex::BCType::int_dir);
            else
                set_bcrec_hi(dir, amrex::BCType::int_dir);
            break;

        case BC::pressure_inflow:
        case BC::pressure_outflow:
        case BC::zero_gradient:
            if (side == amrex::Orientation::low)
                set_bcrec_lo(dir, amrex::BCType::foextrap);
            else
                set_bcrec_hi(dir, amrex::BCType::foextrap);
            break;

        case BC::mass_inflow:
        case BC::no_slip_wall:
            if (side == amrex::Orientation::low)
                set_bcrec_lo(dir, amrex::BCType::ext_dir);
            else
                set_bcrec_hi(dir, amrex::BCType::ext_dir);
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
    auto& fname = m_field.name();
    auto& bctype = m_field.bc_type();
    auto& bcval = m_field.bc_values();
    const int ndim = m_field.num_comp();
    for (amrex::OrientationIter oit; oit; ++oit) {
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
    auto& fname = m_field.name();
    auto& bctype = m_field.bc_type();
    auto& bcval = m_field.bc_values();
    const int ndim = m_field.num_comp();
    for (amrex::OrientationIter oit; oit; ++oit) {
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
    auto& ibctype = m_field.bc_type();
    for (amrex::OrientationIter oit; oit; ++oit) {
        auto ori = oit();
        const auto side = ori.faceDir();
        const auto bct = ibctype[ori];
        const int dir = ori.coordDir();

        switch (bct) {
        case BC::periodic:
            if (side == amrex::Orientation::low)
                set_bcrec_lo(dir, amrex::BCType::int_dir);
            else
                set_bcrec_hi(dir, amrex::BCType::int_dir);
            break;

        default:
            if (side == amrex::Orientation::low)
                set_bcrec_lo(dir, amrex::BCType::foextrap);
            else
                set_bcrec_hi(dir, amrex::BCType::foextrap);
            break;
        }
    }
}

} // namespace amr_wind
