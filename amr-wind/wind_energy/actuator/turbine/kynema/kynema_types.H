#ifndef KYNEMA_TYPES_H
#define KYNEMA_TYPES_H

#include "amr-wind/wind_energy/actuator/turbine/kynema/kynema_wrapper.H"
#include "amr-wind/wind_energy/actuator/turbine/external/external_base_types.H"

#include <string>

namespace ext_turb {

/** Representation of a turbine for exchanging data with Kynema
 */
struct KynemaTurbine : public ext_turb::ExternalTurbine
{
    /** Variables exclusive to Kynema turbine type
     */
    std::unique_ptr<kynema::interfaces::TurbineInterface> interface;

    // Need buffers to make data contiguous to fit into coupling
    void allocate_buffers()
    {
        buffer_cf.resize(length_force(0));
        buffer_pv.resize(
            length_fluid_velocity(0) + length_fluid_velocity(1) +
            length_fluid_velocity(2));
        buffer_sv.resize(
            length_fluid_velocity(0) + length_fluid_velocity(1) +
            length_fluid_velocity(2));
        buffer_fv.resize(
            length_fluid_velocity(0) + length_fluid_velocity(1) +
            length_fluid_velocity(2));
        buffer_fv_2idx.resize(length_fluid_velocity(0));

        buffer_f.resize(length_force(0) + length_force(1) + length_force(2));
        buffer_pf.resize(
            length_position_at_force(0) + length_position_at_force(1) +
            length_position_at_force(2));
        buffer_ori.resize(length_orientation());

        ptr_cf = &buffer_cf[0];
        ptr_pv = &buffer_pv[0];
        ptr_sv = &buffer_sv[0];
        ptr_fv = &buffer_fv[0];
        ptr_f = &buffer_f[0];
        ptr_pf = &buffer_pf[0];
        ptr_ori = &buffer_ori[0];
    }

    void populate_buffers()
    {
        int idx_f = 0;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            // Fluid velocity (fv) is not populated here, it comes from the
            // CFD solver. Need to add step to pass it from buffer to Kynema
            for (int n = 0; n < length_force(0); ++n) {
                if (dir == 0) {
                    buffer_cf[n] = (float)interface->Aerodynamics()
                                       .bodies[body_f_index(n)]
                                       .chord[node_f_index(n)];
                }
                buffer_f[idx_f] = (float)interface->Aerodynamics()
                                      .bodies[body_f_index(n)]
                                      .loads[node_f_index(n)][dir];
                buffer_pf[idx_f] = (float)interface->Aerodynamics()
                                       .bodies[body_f_index(n)]
                                       .x_motion[node_f_index(n)][dir];
                buffer_pv[idx_f] = buffer_pf[idx_f];
                buffer_sv[idx_f] = (float)interface->Aerodynamics()
                                       .bodies[body_f_index(n)]
                                       .v_motion[node_f_index(n)][dir];
                ++idx_f;
            }
        }
        // Orientation stores a 3x3 tensor (9 elements) for every point
        std::array<float, 4> q;
        for (int n = 0; n < length_force(0); ++n) {
            // convert quaternions to rotation matrix
            for (int iq = 0; iq < 4; ++iq) {
                q[iq] = (float)interface->Aerodynamics()
                            .bodies[body_f_index(n)]
                            .qqr_motion_map[node_f_index(n)][iq];
            }
            const int idx = n * 9;
            buffer_ori[idx + 0] =
                q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3];
            buffer_ori[idx + 1] = 2. * (q[1] * q[2] - q[0] * q[3]);
            buffer_ori[idx + 2] = 2. * (q[1] * q[3] + q[0] * q[2]);
            buffer_ori[idx + 3] = 2. * (q[1] * q[2] + q[0] * q[3]);
            buffer_ori[idx + 4] =
                (q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3]);
            buffer_ori[idx + 5] = 2. * (q[2] * q[3] - q[0] * q[1]);
            buffer_ori[idx + 6] = 2. * (q[1] * q[3] - q[0] * q[2]);
            buffer_ori[idx + 7] = 2. * (q[2] * q[3] + q[0] * q[1]);
            buffer_ori[idx + 8] =
                (q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3]);
        }
    }

    void pass_fluid_velocity_directly()
    {
        int idx_v = 0;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            for (int n = 0; n < length_fluid_velocity(0); ++n) {
                interface->Aerodynamics()
                    .bodies[body_f_index(n)]
                    .v_inflow[node_f_index(n)][dir] = buffer_fv[idx_v];
                ++idx_v;
            }
        }
    }

    /** Helper functions for navigating indices */
    int body_f_index(int id) const
    {
        int b_id = 0;
        int id_diff = id;
        for (int b = 0; b < (int)interface->Aerodynamics().bodies.size(); ++b) {
            id_diff -= static_cast<int>(
                interface->Aerodynamics().bodies[b].loads.size());
            if (id_diff < 0) {
                b_id = b;
                break;
            }
        }
        return b_id;
    }

    int node_f_index(int id) const
    {
        int n_id = 0;
        int id_diff = id;
        for (int b = 0; b < (int)interface->Aerodynamics().bodies.size(); ++b) {
            id_diff -= (int)interface->Aerodynamics().bodies[b].loads.size();
            if (id_diff < 0) {
                n_id = id_diff +
                       (int)interface->Aerodynamics().bodies[b].loads.size();
                break;
            }
        }
        return n_id;
    }

    /** Data access functions needed for every turbine type */

    float* position_at_vel(int dir) const override
    {
        return ptr_pv + dir * length_fluid_velocity(0);
    }

    float* solid_velocity(int dir) const override
    {
        return ptr_sv + dir * length_fluid_velocity(0);
    }

    float* fluid_velocity(int dir) const override
    {
        return ptr_fv + dir * length_fluid_velocity(0);
    }

    float* force(int dir) const override
    {
        return ptr_f + dir * length_force(0);
    }

    float* position_at_force(int dir) const override
    {
        return ptr_pf + dir * length_position_at_force(0);
    }

    float* chord_at_force() const override { return ptr_cf; }

    float* orientation() const override { return ptr_ori; }

    int length_fluid_velocity(int dir) const override
    {
        return length_force(dir);
    }

    int length_force(int dir) const override
    {
        unsigned int n_sect{0};
        for (unsigned int n = 0; n < interface->Aerodynamics().bodies.size();
             ++n) {
            n_sect += interface->Aerodynamics().bodies[n].loads.size();
        }
        return static_cast<int>(n_sect);
    }

    int length_position_at_force(int dir) const override
    {
        return length_force(dir);
    }

    int length_orientation() const override
    {
        return length_force(0) * AMREX_SPACEDIM * AMREX_SPACEDIM;
    }

    double fluid_density{1.};

    std::vector<float> buffer_cf;
    std::vector<float> buffer_pv;
    std::vector<float> buffer_sv;
    std::vector<float> buffer_fv;
    std::vector<float> buffer_f;
    std::vector<float> buffer_pf;
    std::vector<float> buffer_ori;
    float* ptr_cf{nullptr};
    float* ptr_pv{nullptr};
    float* ptr_sv{nullptr};
    float* ptr_fv{nullptr};
    float* ptr_f{nullptr};
    float* ptr_pf{nullptr};
    float* ptr_ori{nullptr};

    std::vector<std::array<double, 3>> buffer_fv_2idx;
};

struct KynemaSolverData
{
    std::string m_output_dir{"kynema_velocity_data"};
};

} // namespace ext_turb

#endif /* KYNEMA_TYPES_H */
