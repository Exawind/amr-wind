#ifndef INCFLO_ENUMS_H
#define INCFLO_ENUMS_H

enum struct DiffusionType { Invalid, Explicit, Crank_Nicolson, Implicit };

enum struct BC {
    pressure_outflow,
    mass_inflow,
    mass_inflow_outflow,
    zero_gradient,
    no_slip_wall,
    slip_wall,
    wall_model,
    fixed_gradient,
    periodic,
    wave_generation,
    symmetric_wall,
    undefined
};

#endif /* INCFLO_ENUMS_H */
