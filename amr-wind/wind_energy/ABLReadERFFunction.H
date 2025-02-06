#ifndef ABLREADERFFUNCTION_H
#define ABLREADERFFUNCTION_H
#include <functional>

class MultiBlockContainer;
namespace amr_wind {
class InletData;
}

using ReadERFFunction = std::function<void(
    const amrex::Real time,
    amrex::Vector<amrex::Real>&,
    amr_wind::InletData&,
    const amrex::Vector<amr_wind::Field*>&,
    MultiBlockContainer*)>;

#endif /* ABLREADERFFUNCTION_H */
