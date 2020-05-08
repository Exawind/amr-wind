#include "incflo.H"
#include "FieldFillPatchOps.H"
#include "FieldBCOps.H"
#include "PDE.H"
#include "SchemeTraits.H"

void incflo::declare_fields()
{
    auto& pde_mgr = m_sim.pde_manager();

    // Always register incompressible Navier-Stokes equation
    pde_mgr.register_icns();

    // Register density first so that we can compute its `n+1/2` state before
    // other scalars attempt to use it in their computations.
    if (!m_constant_density) {
        if (pde_mgr.scalar_eqns().size() > 0)
            amrex::Abort(
                "For non-constant density, it must be the first equation "
                "registered for the scalar equations");
        pde_mgr.register_transport_pde("Density");
    }

    // TODO: This should be customized based on Physics
    // pde_mgr.register_transport_pde("Temperature");

    m_sim.create_turbulence_model();
    m_sim.init_physics();
}

void incflo::init_field_bcs ()
{
    using namespace amrex;

    // store temporary array to copy into each field below
    amrex::GpuArray<BC, AMREX_SPACEDIM*2> bc_temp;

    auto f = [&] (std::string const& bcid, Orientation ori)
    {
        ParmParse pp(bcid);
        std::string bc_type_in = "null";
        pp.query("type", bc_type_in);
        std::string bc_type = amrex::toLower(bc_type_in);

        if (bc_type == "pressure_inflow" or bc_type == "pi")
        {
            bc_temp[ori] = BC::pressure_inflow;
        }
        else if (bc_type == "pressure_outflow" or bc_type == "po")
        {
            bc_temp[ori] = BC::pressure_outflow;
        }
        else if (bc_type == "mass_inflow" or bc_type == "mi")
        {
            bc_temp[ori] = BC::mass_inflow;
        }
        else if (bc_type == "no_slip_wall" or bc_type == "nsw")
        {
            bc_temp[ori] = BC::no_slip_wall;
        }
        else if (bc_type == "slip_wall" or bc_type == "sw")
        {
            bc_temp[ori] = BC::slip_wall;
        }
        else if (bc_type == "wall_model" or bc_type == "wm")
        {
            bc_temp[ori] = BC::wall_model;
        }
        else
        {
            bc_temp[ori] = BC::undefined;
        }

        if (geom[0].isPeriodic(ori.coordDir())) {
            if (bc_temp[ori] == BC::undefined) {
                bc_temp[ori] = BC::periodic;
            } else {
                amrex::Abort("Wrong BC type for periodic boundary");
            }
        }
    };

    f("xlo", Orientation(Direction::x,Orientation::low));
    f("xhi", Orientation(Direction::x,Orientation::high));
    f("ylo", Orientation(Direction::y,Orientation::low));
    f("yhi", Orientation(Direction::y,Orientation::high));
    f("zlo", Orientation(Direction::z,Orientation::low));
    f("zhi", Orientation(Direction::z,Orientation::high));

    icns().init_bcs(bc_temp);
    for (auto& eqn: scalar_eqns()) {
        eqn->init_bcs(bc_temp);
    }
}
