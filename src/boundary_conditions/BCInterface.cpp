#include "BCInterface.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace {
amrex::Vector<std::string> bcnames = {"xlo", "ylo", "zlo", "xhi", "yhi", "zhi"};
}

BCIface::BCIface(Field& field, const IncfloBC& ibc_type)
    : m_field(field), m_ibc_type(ibc_type)
{}

inline void BCIface::set_incflo_bc()
{
    auto& ibctype = m_field.bc_type();
    for (int i = 0; i < AMREX_SPACEDIM * 2; ++i) ibctype[i] = m_ibc_type[i];
}

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
    set_incflo_bc();
    set_default_value(value);
    set_bcrec();
    read_values();
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

void BCVelocity::set_bcrec()
{
    auto& bcrec = m_field.bcrec();

    for (amrex::OrientationIter oit; oit; ++oit) {
        auto ori = oit();
        const auto side = ori.faceDir();
        const auto bct = m_ibc_type[ori];
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
            break;
        }
    }
}

void BCScalar::set_bcrec()
{
    for (amrex::OrientationIter oit; oit; ++oit) {
        auto ori = oit();
        const auto side = ori.faceDir();
        const auto bct = m_ibc_type[ori];
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
        case BC::no_slip_wall:
            if (side == amrex::Orientation::low)
                set_bcrec_lo(dir, amrex::BCType::foextrap);
            else
                set_bcrec_hi(dir, amrex::BCType::foextrap);
            break;

        case BC::mass_inflow:
            if (side == amrex::Orientation::low)
                set_bcrec_lo(dir, amrex::BCType::ext_dir);
            else
                set_bcrec_hi(dir, amrex::BCType::ext_dir);
            break;

        case BC::slip_wall:
        case BC::wall_model:
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
        case BC::slip_wall:
        case BC::wall_model:
            pp.queryarr(fname.c_str(), bcval[ori], 0, ndim);
            break;

        default:
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
    for (amrex::OrientationIter oit; oit; ++oit) {
        auto ori = oit();
        const auto side = ori.faceDir();
        const auto bct = m_ibc_type[ori];
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
