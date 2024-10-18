#ifndef ABLREADERFFCN_H
#define ABLREADERFFCN_H
#include <functional>

class MultiBlockContainer;
namespace amr_wind {
class InletData;
}

using ReadERFFunction = std::function<void(
    const amr_wind::SimTime&,
    amrex::Vector<amrex::Real>&,
    amr_wind::InletData&,
    const amrex::Vector<amr_wind::Field*>&,
    MultiBlockContainer*)>;

#endif /* ABLREADERFFCN_H */
