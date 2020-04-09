#include "incflo.H"
#include "FieldFillPatchOps.H"
#include "FieldBCOps.H"
#include "prob_bc.H"

void incflo::declare_fields()
{
    const int nstates = 2;
    const int ng = nghost_state();

    auto& vel = m_repo.declare_cc_field("velocity", AMREX_SPACEDIM, ng, nstates);
    auto& rho = m_repo.declare_cc_field("density", 1, ng, nstates);
    auto& trac = m_repo.declare_cc_field("tracer", m_ntrac, ng, nstates);
    auto& gp = m_repo.declare_cc_field("gp", AMREX_SPACEDIM, 0, 1);
    auto& p = m_repo.declare_nd_field("p", 1, 0, 1);

    auto& vel_for = m_repo.declare_cc_field("velocity_forces", AMREX_SPACEDIM, nghost_force(), 1);
    auto& tra_for = m_repo.declare_cc_field("tracer_forces", m_ntrac, nghost_force(), 1);

    m_repo.declare_cc_field("viscosity", 1, 1, 1);
    m_repo.declare_cc_field("tracer_viscosity", m_ntrac, 1, 1);

    m_repo.declare_cc_field("conv_velocity", AMREX_SPACEDIM, 0, nstates);
    m_repo.declare_cc_field("conv_density", 1, 0, nstates);
    m_repo.declare_cc_field("conv_tracer", m_ntrac, 0, nstates);

    m_repo.declare_cc_field("divtau", AMREX_SPACEDIM, 0, nstates);
    m_repo.declare_cc_field("laps", m_ntrac, 0, nstates);

    m_repo.declare_face_normal_field({"u_mac","v_mac","w_mac"}, 1, nghost_mac(), 1);

    vel.register_fill_patch_op<amr_wind::FieldFillPatchOps<amr_wind::FieldBCDirichlet>>(
        *this, m_time, m_probtype);
    rho.register_fill_patch_op<amr_wind::FieldFillPatchOps<amr_wind::FieldBCDirichlet>>(
        *this, m_time, m_probtype);
    trac.register_fill_patch_op<amr_wind::FieldFillPatchOps<amr_wind::FieldBCDirichlet>>(
        *this, m_time, m_probtype);
    gp.register_fill_patch_op<amr_wind::FieldFillPatchOps<amr_wind::FieldBCNoOp>>(
        *this, m_time, m_probtype);

    p.register_fill_patch_op<amr_wind::FieldFillConstScalar>(0.0);

    vel_for.register_fill_patch_op<amr_wind::FieldFillPatchOps<amr_wind::FieldBCNoOp>>(
        *this, m_time, m_probtype, amr_wind::FieldInterpolator::PiecewiseConstant);
    tra_for.register_fill_patch_op<amr_wind::FieldFillPatchOps<amr_wind::FieldBCNoOp>>(
        *this, m_time, m_probtype, amr_wind::FieldInterpolator::PiecewiseConstant);

    // Inform field repo which fields need fillpatch operations on regrid
    vel.fillpatch_on_regrid() = true;
    rho.fillpatch_on_regrid() = true;
    trac.fillpatch_on_regrid() = true;
    gp.fillpatch_on_regrid() = true;

}

void incflo::init_field_bcs ()
{
    using namespace amrex;
    auto& velocity = m_repo.get_field("velocity");
    auto& density = m_repo.get_field("density");
    auto& tracer = m_repo.get_field("tracer");
    auto& vel_for = m_repo.get_field("velocity_forces");
    auto& tra_for = m_repo.get_field("tracer_forces");

    auto& bc_velocity = velocity.bc_values();
    auto& bcrec_velocity = velocity.bcrec();
    auto& bc_density = density.bc_values();
    auto& bcrec_density = density.bcrec();
    auto& bc_tracer = tracer.bc_values();
    auto& bcrec_tracer = tracer.bcrec();
    auto& bcrec_vel_for = vel_for.bcrec();
    auto& bcrec_tra_for = tra_for.bcrec();

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

        if (bc_type == "mass_inflow" or bc_type == "mi")
        {
            std::vector<Real> v;
            if (pp.queryarr("velocity", v, 0, AMREX_SPACEDIM)) {
               bc_velocity[ori] = {v[0],v[1],v[2]};
            }

            pp.queryarr("density", bc_density[ori], 0, 1);

            pp.queryarr("tracer", bc_tracer[ori], 0, m_ntrac);
        }
        else if (bc_type == "no_slip_wall" or bc_type == "nsw")
        {
            std::vector<Real> v;
            if (pp.queryarr("velocity", v, 0, AMREX_SPACEDIM)) {
                v[ori.coordDir()] = 0.0;
                bc_velocity[ori] = {v[0],v[1],v[2]};
            }
        }
        else if (bc_type == "slip_wall" or bc_type == "sw")
        {
            pp.queryarr("tracer", bc_tracer[ori], 0, m_ntrac);

        }
        else if (bc_type == "wall_model" or bc_type == "wm")
        {
            pp.queryarr("tracer", bc_tracer[ori], 0, m_ntrac);
        }
    };

    f("xlo", Orientation(Direction::x,Orientation::low));
    f("xhi", Orientation(Direction::x,Orientation::high));
    f("ylo", Orientation(Direction::y,Orientation::low));
    f("yhi", Orientation(Direction::y,Orientation::high));
    f("zlo", Orientation(Direction::z,Orientation::low));
    f("zhi", Orientation(Direction::z,Orientation::high));

    {
        for (OrientationIter oit; oit; ++oit) {
            Orientation ori = oit();
            int dir = ori.coordDir();
            Orientation::Side side = ori.faceDir();
            auto const bct = m_bc_type[ori];
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
            auto const bct = m_bc_type[ori];
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
            auto const bct = m_bc_type[ori];
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
            auto const bct = m_bc_type[ori];
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
            auto const bct = m_bc_type[ori];
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
