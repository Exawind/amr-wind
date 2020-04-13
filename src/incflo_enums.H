#ifndef INCFLO_ENUMS_H
#define INCFLO_ENUMS_H

enum struct FluidModel {
    Newtonian,
    powerlaw,
    Bingham,
    HerschelBulkley,
    deSouzaMendesDutra,
    SmagorinskyLillySGS
};

enum struct DiffusionType { Invalid, Explicit, Crank_Nicolson, Implicit };

enum struct BC {
    pressure_inflow,
    pressure_outflow,
    mass_inflow,
    no_slip_wall,
    slip_wall,
    wall_model,
    periodic,
    undefined
};

#endif /* INCFLO_ENUMS_H */
