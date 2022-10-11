#ifndef TRANSPORTMODEL_H
#define TRANSPORTMODEL_H

#include "amr-wind/core/Factory.H"
#include "amr-wind/core/FieldRepo.H"

namespace amr_wind::transport {

/**
 *  \defgroup transport Transport models
 *
 *  AMR-Wind representation of fluid transport models.
 *
 *  For most land-based wind-energy applications, a constant transport
 *  (ConstTransport) model is sufficient with uniform properties. However, for
 *  offshore applications with multiphase models, the properties can change at
 *  different positions in the domain. Transport models provide a way to
 *  separate out this behavior from the base turbulence model implementation.
 *  The turbulence models are templated on transport model to provide different
 *  behavior depending on fluid transport model chosen by the user. For more
 *  details see amr_wind::turbulence::TurbModel implementation.
 *
 */

/** Abstract representation of a transport model
 *  \ingroup transport
 */
class TransportModel
{
public:
    virtual ~TransportModel() = default;

    //! Dynamic laminar viscosity (kg/m/s)
    virtual std::unique_ptr<ScratchField> mu() = 0;

    //! Thermal diffusivity
    virtual std::unique_ptr<ScratchField> alpha() = 0;

    //! Scalar diffusivity based on Schmidt number
    virtual std::unique_ptr<ScratchField>
    scalar_diffusivity(const std::string& scalar_name) = 0;
};
} // namespace amr_wind::transport

#endif /* TRANSPORTMODEL_H */
