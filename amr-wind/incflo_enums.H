#ifndef INCFLO_ENUMS_H
#define INCFLO_ENUMS_H

enum struct DiffusionType { Invalid, Explicit, Crank_Nicolson, Implicit };

enum struct BC {
    pressure_inflow,
    pressure_outflow,
    mass_inflow,
    zero_gradient,
    no_slip_wall,
    slip_wall,
    wall_model,
    fixed_gradient,
    periodic,
    wave_generation,
    undefined
};

#endif /* INCFLO_ENUMS_H */
