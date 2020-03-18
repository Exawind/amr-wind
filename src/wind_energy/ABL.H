#ifndef ABL_H
#define ABL_H

#include <memory>

#include "ABLFieldInit.H"

namespace amr_wind {

/** Atmospheric Boundary Layer physics
 */
class ABL
{
public:
    ABL();

    void add_abl_source_terms();

    const ABLFieldInit& abl_initializer() const { return *m_field_init; }

private:
    //! ABL field initializer instance
    std::unique_ptr<ABLFieldInit> m_field_init;

    //! True if Boussinesq buoyancy term is active
    bool m_has_boussinesq{true};

    //! True if Coriolis forcing term is active
    bool m_has_coriolis{true};

    //! True if driving pressure gradient term is activate
    bool m_has_driving_dpdx{true};
};

}

#endif /* ABL_H */
