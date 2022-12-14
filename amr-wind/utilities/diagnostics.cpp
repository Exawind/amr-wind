#include "amr-wind/incflo.H"

using namespace amrex;

//
// Print maximum values (useful for tracking evolution)
//
void incflo::PrintMaxValues(const std::string& header)
{
    BL_PROFILE("amr-wind::incflo::PrintMaxValues");

    amrex::Print() << "\nL-inf norm summary: " << header << std::endl
                   << "........................................................"
                      "......................";

    for (int lev = 0; lev <= finest_level; lev++) {
        amrex::Print() << "\nLevel " << lev << std::endl;

        const auto& vel = icns().fields().field;
        amrex::Print() << "  " << std::setw(16) << std::left << vel.name();
        for (int i = 0; i < vel.num_comp(); ++i) {
            amrex::Print() << std::setw(20) << std::right << vel(lev).norm0(i);
        }
        amrex::Print() << std::endl;

        const auto& gradp = grad_p();
        amrex::Print() << "  " << std::setw(16) << std::left << gradp.name();
        for (int i = 0; i < gradp.num_comp(); ++i) {
            amrex::Print() << std::setw(20) << std::right
                           << gradp(lev).norm0(i);
        }
        amrex::Print() << std::endl;

        for (auto& eqn : scalar_eqns()) {
            auto& field = eqn->fields().field;
            amrex::Print() << "  " << std::setw(16) << std::left
                           << field.name();
            for (int i = 0; i < field.num_comp(); ++i) {
                amrex::Print()
                    << std::setw(20) << std::right << field(lev).norm0(i);
            }
            amrex::Print() << std::endl;
        }
    }
    amrex::Print() << "........................................................"
                      "......................"
                   << std::endl
                   << std::endl;
}

//
// Print the maximum values of the velocity components and velocity divergence
//
void incflo::PrintMaxVel(int lev) const
{
    BL_PROFILE("amr-wind::incflo::PrintMaxVel");
    amrex::Print() << "max(abs(u/v/w))  = " << velocity()(lev).norm0(0) << "  "
                   << velocity()(lev).norm0(1) << "  "
                   << velocity()(lev).norm0(2) << "  " << std::endl;
}

//
// Print the maximum values of the pressure gradient components and pressure
//
void incflo::PrintMaxGp(int lev) const
{
    BL_PROFILE("amr-wind::incflo::PrintMaxGp");
    amrex::Print() << "max(abs(gpx/gpy/gpz/p))  = " << grad_p()(lev).norm0(0)
                   << "  " << grad_p()(lev).norm0(1) << "  "
                   << grad_p()(lev).norm0(2) << "  " << pressure()(lev).norm0(0)
                   << "  " << std::endl;
}

void incflo::CheckForNans(int lev) const
{
    BL_PROFILE("amr-wind::incflo::CheckForNans");
    bool ro_has_nans = density()(lev).contains_nan(false);
    bool ug_has_nans = velocity()(lev).contains_nan(false);
    bool vg_has_nans = velocity()(lev).contains_nan(true);
    bool wg_has_nans = velocity()(lev).contains_nan(true);
    bool pg_has_nans = pressure()(lev).contains_nan(false);

    if (ro_has_nans) {
        amrex::Print() << "WARNING: ro contains NaNs!!!";
    }

    if (ug_has_nans) {
        amrex::Print() << "WARNING: u contains NaNs!!!";
    }

    if (vg_has_nans) {
        amrex::Print() << "WARNING: v contains NaNs!!!";
    }

    if (wg_has_nans) {
        amrex::Print() << "WARNING: w contains NaNs!!!";
    }

    if (pg_has_nans) {
        amrex::Print() << "WARNING: p contains NaNs!!!";
    }
}
