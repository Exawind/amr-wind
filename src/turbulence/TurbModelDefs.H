#ifndef TURBMODELDEFS_H
#define TURBMODELDEFS_H

#include "TurbulenceModel.H"
#include "ConstTransport.H"

namespace amr_wind {

#define INSTANTIATE_TURBULENCE_MODEL(Model) \
    template struct ::amr_wind::turbulence::TurbulenceModel::Register<::amr_wind::turbulence::Model<transport::ConstTransport>>

}

#endif /* TURBMODELDEFS_H */
