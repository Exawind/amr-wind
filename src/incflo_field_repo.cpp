#include "incflo.H"
#include "FieldFillPatchOps.H"
#include "FieldBCOps.H"
#include "prob_bc.H"
#include "PDE.H"
#include "SchemeTraits.H"

void incflo::declare_fields()
{
    const std::string scheme = m_use_godunov
                                   ? amr_wind::fvm::Godunov::scheme_name()
                                   : amr_wind::fvm::MOL::scheme_name();
    m_icns = amr_wind::pde::PDEBase::create(
        "ICNS-" + scheme, m_time, m_repo, m_probtype);
    m_scalar_eqns.emplace_back(amr_wind::pde::PDEBase::create(
        "Temperature-" + scheme, m_time, m_repo, m_probtype));
    m_scalar_eqns.emplace_back(amr_wind::pde::PDEBase::create(
        "Density-" + scheme, m_time, m_repo, m_probtype));
}

void incflo::init_field_bcs ()
{
    using namespace amrex;
    auto& velocity = m_repo.get_field("velocity");
    auto& density = m_repo.get_field("density");
    auto& tracer = m_repo.get_field("temperature");
    auto& vel_for = m_repo.get_field("velocity_src_term");
    auto& tra_for = m_repo.get_field("temperature_src_term");

    auto& bc_velocity = velocity.bc_values();
    auto& bcrec_velocity = velocity.bcrec();
    auto& bc_density = density.bc_values();
    auto& bcrec_density = density.bcrec();
    auto& bc_tracer = tracer.bc_values();
    auto& bcrec_tracer = tracer.bcrec();
    auto& bcrec_vel_for = vel_for.bcrec();
    auto& bcrec_tra_for = tra_for.bcrec();
    auto& bc_pressure = pressure().bc_values();

    // store temporary array to copy into each field below
    amrex::GpuArray<BC, AMREX_SPACEDIM*2> bc_temp;

    auto f = [&] (std::string const& bcid, Orientation ori)
    {
        bc_density[ori][0] = 1.0;
        bc_velocity[ori][0] = 0.0; // default
        bc_velocity[ori][1] = 0.0;
        bc_velocity[ori][2] = 0.0;

        ParmParse pp(bcid);
        std::string bc_type_in = "null";
        pp.query("type", bc_type_in);
        std::string bc_type = amrex::toLower(bc_type_in);

        if (bc_type == "pressure_inflow" or bc_type == "pi")
        {
            bc_temp[ori] = BC::pressure_inflow;
            pp.get("pressure", bc_pressure[ori][0]);
        }
        else if (bc_type == "pressure_outflow" or bc_type == "po")
        {
            bc_temp[ori] = BC::pressure_outflow;
            pp.get("pressure", bc_pressure[ori][0]);
        }
        else if (bc_type == "mass_inflow" or bc_type == "mi")
        {
            bc_temp[ori] = BC::mass_inflow;

            std::vector<Real> v;
            if (pp.queryarr("velocity", v, 0, AMREX_SPACEDIM)) {
               bc_velocity[ori] = {v[0],v[1],v[2]};
            }

            pp.queryarr("density", bc_density[ori], 0, 1);

            pp.queryarr("tracer", bc_tracer[ori], 0, m_ntrac);
        }
        else if (bc_type == "no_slip_wall" or bc_type == "nsw")
        {
            bc_temp[ori] = BC::no_slip_wall;

            std::vector<Real> v;
            if (pp.queryarr("velocity", v, 0, AMREX_SPACEDIM)) {
                v[ori.coordDir()] = 0.0;
                bc_velocity[ori] = {v[0],v[1],v[2]};
            }
        }
        else if (bc_type == "slip_wall" or bc_type == "sw")
        {
            bc_temp[ori] = BC::slip_wall;
            pp.queryarr("tracer", bc_tracer[ori], 0, m_ntrac);
        }
        else if (bc_type == "wall_model" or bc_type == "wm")
        {
            bc_temp[ori] = BC::wall_model;
            m_wall_model_flag = true;
            pp.queryarr("tracer", bc_tracer[ori], 0, m_ntrac);
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

    for (int i=0; i < AMREX_SPACEDIM*2; ++i) {
        velocity.bc_type()[i] = bc_temp[i];
        tracer.bc_type()[i] = bc_temp[i];
        density.bc_type()[i] = bc_temp[i];
        pressure().bc_type()[i] = bc_temp[i];
    }

    {
        for (OrientationIter oit; oit; ++oit) {
            Orientation ori = oit();
            int dir = ori.coordDir();
            Orientation::Side side = ori.faceDir();
            auto const bct = bc_temp[ori];
            if (bct == BC::pressure_inflow or
                bct == BC::pressure_outflow)
            {
                if (side == Orientation::low) {
                    bcrec_velocity[0].setLo(dir, BCType::foextrap);
                    bcrec_velocity[1].setLo(dir, BCType::foextrap);
                    bcrec_velocity[2].setLo(dir, BCType::foextrap);
                } else {
                    bcrec_velocity[0].setHi(dir, BCType::foextrap);
                    bcrec_velocity[1].setHi(dir, BCType::foextrap);
                    bcrec_velocity[2].setHi(dir, BCType::foextrap);
                }
            }
            else if (bct == BC::mass_inflow or bct == BC::no_slip_wall)
            {
                if (side == Orientation::low) {
                    bcrec_velocity[0].setLo(dir, BCType::ext_dir);
                    bcrec_velocity[1].setLo(dir, BCType::ext_dir);
                    bcrec_velocity[2].setLo(dir, BCType::ext_dir);
                } else {
                    bcrec_velocity[0].setHi(dir, BCType::ext_dir);
                    bcrec_velocity[1].setHi(dir, BCType::ext_dir);
                    bcrec_velocity[2].setHi(dir, BCType::ext_dir);
                }
            }
            else if (bct == BC::slip_wall or bct == BC::wall_model) //fixme not sure what to do yet for wall model... first order or high order?
            {
                if (side == Orientation::low) {
                    // Tangential directions have hoextrap
                    bcrec_velocity[0].setLo(dir, BCType::hoextrap);
                    bcrec_velocity[1].setLo(dir, BCType::hoextrap);
                    bcrec_velocity[2].setLo(dir, BCType::hoextrap);

                    // Only normal direction has ext_dir
                    bcrec_velocity[dir].setLo(dir, BCType::ext_dir);

                } else {
                    // Tangential directions have hoextrap
                    bcrec_velocity[0].setHi(dir, BCType::hoextrap);
                    bcrec_velocity[1].setHi(dir, BCType::hoextrap);
                    bcrec_velocity[2].setHi(dir, BCType::hoextrap);

                    // Only normal direction has ext_dir
                    bcrec_velocity[dir].setHi(dir, BCType::ext_dir);
                }
            }
            else if (bct == BC::periodic)
            {
                if (side == Orientation::low) {
                    bcrec_velocity[0].setLo(dir, BCType::int_dir);
                    bcrec_velocity[1].setLo(dir, BCType::int_dir);
                    bcrec_velocity[2].setLo(dir, BCType::int_dir);
                } else {
                    bcrec_velocity[0].setHi(dir, BCType::int_dir);
                    bcrec_velocity[1].setHi(dir, BCType::int_dir);
                    bcrec_velocity[2].setHi(dir, BCType::int_dir);
                }
            }
        }
    }

    {
        for (OrientationIter oit; oit; ++oit) {
            Orientation ori = oit();
            int dir = ori.coordDir();
            Orientation::Side side = ori.faceDir();
            auto const bct = bc_temp[ori];
            if (bct == BC::pressure_inflow  or
                bct == BC::pressure_outflow or
                bct == BC::no_slip_wall)
            {
                if (side == Orientation::low) {
                    bcrec_density[0].setLo(dir, BCType::foextrap);
                } else {
                    bcrec_density[0].setHi(dir, BCType::foextrap);
                }
            }
            else if (bct == BC::slip_wall or
                     bct == BC::wall_model)
            {
                if (side == Orientation::low) {
                    bcrec_density[0].setLo(dir, BCType::hoextrap);
                } else {
                    bcrec_density[0].setHi(dir, BCType::hoextrap);
                }
            }
            else if (bct == BC::mass_inflow)
            {
                if (side == Orientation::low) {
                    bcrec_density[0].setLo(dir, BCType::ext_dir);
                } else {
                    bcrec_density[0].setHi(dir, BCType::ext_dir);
                }
            }
            else if (bct == BC::periodic)
            {
                if (side == Orientation::low) {
                    bcrec_density[0].setLo(dir, BCType::int_dir);
                } else {
                    bcrec_density[0].setHi(dir, BCType::int_dir);
                }
            }
        }
    }

    if (m_ntrac > 0)
    {
        for (OrientationIter oit; oit; ++oit) {
            Orientation ori = oit();
            int dir = ori.coordDir();
            Orientation::Side side = ori.faceDir();
            auto const bct = bc_temp[ori];
            if (bct == BC::pressure_inflow  or
                bct == BC::pressure_outflow or
                bct == BC::no_slip_wall)
            {
                if (side == Orientation::low) {
                    for (auto& b : bcrec_tracer) b.setLo(dir, BCType::foextrap);
                } else {
                    for (auto& b : bcrec_tracer) b.setHi(dir, BCType::foextrap);
                }
            }
            else if (bct == BC::slip_wall   or
                     bct == BC::wall_model)
            {
                if (side == Orientation::low) {
                    for (auto& b : bcrec_tracer) b.setLo(dir, BCType::hoextrap);
                } else {
                    for (auto& b : bcrec_tracer) b.setHi(dir, BCType::hoextrap);
                }
            }
            else if (bct == BC::mass_inflow)
            {
                if (side == Orientation::low) {
                    for (auto& b : bcrec_tracer) b.setLo(dir, BCType::ext_dir);
                } else {
                    for (auto& b : bcrec_tracer) b.setHi(dir, BCType::ext_dir);
                }
            }
            else if (bct == BC::periodic)
            {
                if (side == Orientation::low) {
                    for (auto& b : bcrec_tracer) b.setLo(dir, BCType::int_dir);
                } else {
                    for (auto& b : bcrec_tracer) b.setHi(dir, BCType::int_dir);
                }
            }
        }
    }

    {
        for (OrientationIter oit; oit; ++oit) {
            Orientation ori = oit();
            int dir = ori.coordDir();
            Orientation::Side side = ori.faceDir();
            auto const bct = bc_temp[ori];
            if (bct == BC::periodic)
            {
                if (side == Orientation::low) {
                    for (auto& b : bcrec_vel_for) b.setLo(dir, BCType::int_dir);
                } else {
                    for (auto& b : bcrec_vel_for) b.setHi(dir, BCType::int_dir);
                }
            }
            else
            {
                if (side == Orientation::low) {
                    for (auto& b : bcrec_vel_for) b.setLo(dir, BCType::foextrap);
                } else {
                    for (auto& b : bcrec_vel_for) b.setHi(dir, BCType::foextrap);
                }
            }
        }
    }

    if (m_ntrac > 0)
    {
        for (OrientationIter oit; oit; ++oit) {
            Orientation ori = oit();
            int dir = ori.coordDir();
            Orientation::Side side = ori.faceDir();
            auto const bct = bc_temp[ori];
            if (bct == BC::periodic)
            {
                if (side == Orientation::low) {
                    for (auto& b : bcrec_tra_for) b.setLo(dir, BCType::int_dir);
                } else {
                    for (auto& b : bcrec_tra_for) b.setHi(dir, BCType::int_dir);
                }
            }
            else
            {
                if (side == Orientation::low) {
                    for (auto& b : bcrec_tra_for) b.setLo(dir, BCType::foextrap);
                } else {
                    for (auto& b : bcrec_tra_for) b.setHi(dir, BCType::foextrap);
                }
            }
        }
    }

    velocity.copy_bc_to_device();
    density.copy_bc_to_device();
    tracer.copy_bc_to_device();
    vel_for.copy_bc_to_device();
    tra_for.copy_bc_to_device();

}
