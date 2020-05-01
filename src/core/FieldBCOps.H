#ifndef FIELDBCOPS_H
#define FIELDBCOPS_H

#include "FieldDescTypes.H"

#include "AMReX_Gpu.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_Geometry.H"
#include "AMReX_PhysBCFunct.H"

namespace amr_wind {

class Field;

/** Interface to setup up boundary conditions on a field
 */
class FieldBCIface
{
public:
    virtual ~FieldBCIface() = default;

    virtual void operator()(Field& field, const FieldState rho_state) = 0;
};

struct FieldBCNoOp
{
    int m_probtype;
    int m_ncomp;
    amrex::GpuArray<const amrex::Real*, AMREX_SPACEDIM * 2> m_bcv;

    AMREX_GPU_HOST
    constexpr FieldBCNoOp(
        const int probtype,
        const int ncomp,
        const amrex::GpuArray<const amrex::Real*, AMREX_SPACEDIM * 2>& bcv)
        : m_probtype(probtype), m_ncomp(ncomp), m_bcv(bcv)
    {
    }

    AMREX_GPU_DEVICE
    void operator()(
        const amrex::IntVect& /* iv */,
        amrex::Array4<amrex::Real> const& /* field */,
        const int /* dcomp */,
        const int /* numcomp */,
        amrex::GeometryData const& /* geom */,
        const amrex::Real /* time */,
        const amrex::BCRec* /* bcr */,
        const int /* bcomp */,
        const int /* orig_comp */) const
    {}
};

struct FieldBCDirichlet
{
    int m_probtype;
    int m_ncomp;
    amrex::GpuArray<const amrex::Real*, AMREX_SPACEDIM * 2> m_bcv;

    AMREX_GPU_HOST
    constexpr FieldBCDirichlet(
        const int probtype,
        const int ncomp,
        const amrex::GpuArray<const amrex::Real*, AMREX_SPACEDIM * 2>& bcv)
        : m_probtype(probtype), m_ncomp(ncomp), m_bcv(bcv)
    {
    }

    AMREX_GPU_DEVICE
    void operator()(
        const amrex::IntVect& iv,
        amrex::Array4<amrex::Real> const& field,
        const int /* dcomp */,
        const int /* numcomp */,
        amrex::GeometryData const& geom,
        const amrex::Real /* time */,
        const amrex::BCRec* bcr,
        const int bcomp,
        const int /* orig_comp */) const
    {
        // do something for external Dirichlet (amrex::BCType::ext_dir)
        const int i = iv[0];
        const int j = iv[1];
        const int k = iv[2];

        const amrex::Box& domain_box = geom.Domain();

        for (int n = 0; n < m_ncomp; ++n) {
            const amrex::BCRec& bc = bcr[bcomp + n];

            if (bc.lo(0) == amrex::BCType::ext_dir and
                i < domain_box.smallEnd(0)) {
                field(i, j, k, n) = m_bcv[amrex::Orientation(
                    amrex::Direction::x, amrex::Orientation::low)][n];
            } else if (
                bc.hi(0) == amrex::BCType::ext_dir and
                i > domain_box.bigEnd(0)) {
                field(i, j, k, n) = m_bcv[amrex::Orientation(
                    amrex::Direction::x, amrex::Orientation::high)][n];
            }

            if (bc.lo(1) == amrex::BCType::ext_dir and
                j < domain_box.smallEnd(1)) {
                field(i, j, k, n) = m_bcv[amrex::Orientation(
                    amrex::Direction::y, amrex::Orientation::low)][n];
            } else if (
                bc.hi(1) == amrex::BCType::ext_dir and
                j > domain_box.bigEnd(1)) {
                field(i, j, k, n) = m_bcv[amrex::Orientation(
                    amrex::Direction::y, amrex::Orientation::high)][n];
            }

            if (bc.lo(2) == amrex::BCType::ext_dir and
                k < domain_box.smallEnd(2)) {
                field(i, j, k, n) = m_bcv[amrex::Orientation(
                    amrex::Direction::z, amrex::Orientation::low)][n];
            } else if (
                bc.hi(2) == amrex::BCType::ext_dir and
                k > domain_box.bigEnd(2)) {
                field(i, j, k, n) = m_bcv[amrex::Orientation(
                    amrex::Direction::z, amrex::Orientation::high)][n];
            }
        }
    }
};

} // namespace amr_wind

#endif /* FIELDBCOPS_H */
