#ifndef ABLWALLFUNCTION_H
#define ABLWALLFUNCTION_H

#include "AMReX_REAL.H"

namespace amr_wind {

class ABLWallFunction
{
public:
    ABLWallFunction();

    //! Return the friction velocity given velocity at a height
    amrex::Real utau(const amrex::Real uvel, const amrex::Real zh) const;

private:
    amrex::Real m_kappa{0.41};
    amrex::Real m_z0{0.1};
};

}

#endif /* ABLWALLFUNCTION_H */
