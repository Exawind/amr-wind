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
        ;
        buffer_fv.resize(
            length_fluid_velocity(0) + length_fluid_velocity(1) +
            length_fluid_velocity(2));
        ;
        buffer_f.resize(length_force(0) + length_force(1) + length_force(2));
        ;
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
        int idx = 0;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            for (int n = 0; n < length_force(0); ++n) {
                if (dir == 0) {
                    buffer_cf[n] = (float)interface->Aerodynamics()
                                       .bodies[body_index(n)]
                                       .chord[node_index(n)];
                }
                buffer_pv[idx] = (float)interface->Aerodynamics()
                                     .bodies[body_index(n)]
                                     .node_u[node_index(n)][dir];
                buffer_sv[idx] = (float)interface->Aerodynamics()
                                     .bodies[body_index(n)]
                                     .node_v[node_index(n)][dir];
                // Fluid velocity (fv) is not populated here, it comes from the
                // CFD solver
                buffer_f[idx] = (float)interface->Aerodynamics()
                                    .bodies[body_index(n)]
                                    .node_f[node_index(n)][dir];
                // Is there a distinction between the position of forces and
                // velocities?
                buffer_pf[idx] = buffer_pv[idx];
                ++idx;
            }
        }
        for (int n = 0; n < length_orientation(); ++n) {
            buffer_ori[n] = 0.;
            // put some conversion from node_u to orientation
        }
    }

    /** Helper functions for navigating indices */
    int body_index(int id) const
    {
        int b_id = 0;
        int id_diff = id;
        for (int b = 0; b < (int)interface->Aerodynamics().bodies.size(); ++b) {
            id_diff -= static_cast<int>(
                interface->Aerodynamics().bodies[b].node_ids.size());
            if (id_diff < 0) {
                b_id = b;
                break;
            }
        }
        return b_id;
    }

    int node_index(int id) const
    {
        int n_id = 0;
        int id_diff = id;
        for (int b = 0; b < (int)interface->Aerodynamics().bodies.size(); ++b) {
            id_diff -= (int)interface->Aerodynamics().bodies[b].node_ids.size();
            if (id_diff < 0) {
                n_id = id_diff +
                       (int)interface->Aerodynamics().bodies[b].node_ids.size();
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
        unsigned int n_sect{0};
        for (unsigned int n = 0; n < interface->Aerodynamics().bodies.size();
             ++n) {
            n_sect += interface->Aerodynamics().bodies[n].node_ids.size();
        }
        return static_cast<int>(n_sect);
    }

    int length_force(int dir) const override
    {
        return length_fluid_velocity(dir);
    }

    int length_position_at_force(int dir) const override
    {
        return length_fluid_velocity(dir);
    }

    int length_orientation() const override { return 0; }

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
};

struct KynemaSolverData
{
    std::string m_output_dir{"kynema_velocity_data"};
};

} // namespace ext_turb

#endif /* KYNEMA_TYPES_H */
