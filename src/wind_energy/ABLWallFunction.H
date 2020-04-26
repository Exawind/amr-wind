#ifndef ABLWALLFUNCTION_H
#define ABLWALLFUNCTION_H

#include "CFDSim.H"
#include "PlaneAveraging.H"
#include "FieldBCOps.H"

namespace amr_wind {

class ABLWallFunction
{
public:
    ABLWallFunction(const CFDSim& sim);

    ~ABLWallFunction() = default;

    amrex::Real log_law_height() const { return m_log_law_height; }

    amrex::Real utau() const { return m_utau; }

    const amrex::Array<amrex::Real, AMREX_SPACEDIM>& umean() const
    { return m_umean; }

    void init_log_law_height();

    void update_umean(const PlaneAveraging& pa);

private:
    const amrex::AmrCore& m_mesh;

    amrex::Array<amrex::Real, AMREX_SPACEDIM> m_umean{{0.0, 0.0, 0.0}};

    amrex::Real m_log_law_height{0.0};
    amrex::Real m_utau;

    amrex::Real m_kappa{0.41};
    amrex::Real m_z0{0.1};

    int m_direction{2};
    bool m_use_fch{true};
};

class ABLVelWallFunc : public FieldBCIface
{
public:
    ABLVelWallFunc(Field& velocity, const ABLWallFunction& wall_func);

    void operator()(Field& velocity, const FieldState rho_state) override;

private:
    const ABLWallFunction& m_wall_func;
};

class ABLThetaTopBC : public FieldBCIface
{
public:
    ABLThetaTopBC(Field&) {}

    void operator()(Field& theta, const FieldState rho_state) override;
};

}

#endif /* ABLWALLFUNCTION_H */
