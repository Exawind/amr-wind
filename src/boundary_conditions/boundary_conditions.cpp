#include <AMReX_Vector.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_ParmParse.H>

#include <incflo.H>

using namespace amrex;

void incflo::init_bcs ()
{
    auto f = [this] (std::string const& bcid, Orientation ori)
    {

        ParmParse pp(bcid);
        std::string bc_type_in = "null";
        pp.query("type", bc_type_in);
        std::string bc_type = amrex::toLower(bc_type_in);

        if (bc_type == "pressure_inflow" or bc_type == "pi")
        {
            amrex::Print() << bcid << " set to pressure inflow.\n";

            m_bc_type[ori] = BC::pressure_inflow;

            pp.get("pressure", m_bc_pressure[ori]);
        }
        else if (bc_type == "pressure_outflow" or bc_type == "po")
        {
            amrex::Print() << bcid << " set to pressure outflow.\n";

            m_bc_type[ori] = BC::pressure_outflow;

            pp.get("pressure", m_bc_pressure[ori]);
        }
        else if (bc_type == "mass_inflow" or bc_type == "mi")
        {
            amrex::Print() << bcid << " set to mass inflow.\n";

            m_bc_type[ori] = BC::mass_inflow;
        }
        else if (bc_type == "no_slip_wall" or bc_type == "nsw")
        {
            amrex::Print() << bcid <<" set to no-slip wall.\n";

            m_bc_type[ori] = BC::no_slip_wall;

        }
        else if (bc_type == "slip_wall" or bc_type == "sw")
        {
            amrex::Print() << bcid <<" set to slip wall.\n";

            m_bc_type[ori] = BC::slip_wall;

        }
        else if (bc_type == "wall_model" or bc_type == "wm")
        {
            amrex::Print() << bcid <<" set to wall model.\n";

            m_bc_type[ori] = BC::wall_model;
            m_wall_model_flag = true;

        }
        else
        {
            m_bc_type[ori] = BC::undefined;
        }

        if (geom[0].isPeriodic(ori.coordDir())) {
            if (m_bc_type[ori] == BC::undefined) {
                m_bc_type[ori] = BC::periodic;
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

}

