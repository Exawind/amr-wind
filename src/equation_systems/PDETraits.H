#ifndef PDETRAITS_H
#define PDETRAITS_H

#include <string>
#include "AMReX_MLTensorOp.H"

#include "SourceTerm.H"

namespace amr_wind {
namespace pde {

struct VectorTransport
{
    using MLDiffOp = amrex::MLABecLaplacian;
    static constexpr int ndim = AMREX_SPACEDIM;
};

struct ScalarTransport
{
    using MLDiffOp = amrex::MLABecLaplacian;
    using SrcTerm = SourceTerm;

    static constexpr int ndim = 1;
    static constexpr bool multiply_rho = true;
};

struct ICNS : VectorTransport
{
    using MLDiffOp = amrex::MLTensorOp;
    using SrcTerm = MomentumSource;

    static std::string pde_name() { return "ICNS"; }
    static std::string var_name() { return "velocity"; }
};

struct Density : ScalarTransport
{
    static std::string pde_name() { return "Density"; }
    static std::string var_name() { return "density"; }

    static constexpr bool multiply_rho = false;
};

struct Temperature : ScalarTransport
{
    static std::string pde_name() { return "Temperature"; }
    static std::string var_name() { return "temperature"; }
};

struct TKE : ScalarTransport
{
    static std::string pde_name() { return "TKE"; }
    static std::string var_name() { return "tke"; }
};

struct SDR : ScalarTransport
{
    static std::string pde_name() { return "SDR"; }
    static std::string var_name() { return "sdr"; }
};

} // namespace pde
} // namespace amr_wind

#endif /* PDETRAITS_H */
