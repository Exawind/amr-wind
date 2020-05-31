#ifndef TRANSPORTMODEL_H
#define TRANSPORTMODEL_H

#include "amr-wind/core/Factory.H"
#include "amr-wind/core/FieldRepo.H"

namespace amr_wind {
namespace transport {

class TransportModel
{
public:
    static constexpr bool constant_properties = false;

    virtual ~TransportModel() = default;

    //! Dynamic laminar viscosity (kg/m/s)
    virtual std::unique_ptr<ScratchField> mu() = 0;

    //! Thermal diffusivity
    virtual std::unique_ptr<ScratchField> alpha() = 0;

    //! Scalar diffusivity based on Schmidt number
    virtual std::unique_ptr<ScratchField>
    scalar_diffusivity(const std::string& scalar_name) = 0;
};
}
}

#endif /* TRANSPORTMODEL_H */
