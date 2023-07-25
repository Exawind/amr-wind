#ifndef CHECKPOINTFILEUTIL_H_
#define CHECKPOINTFILEUTIL_H_
#include <AMReX_Config.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>
#include <string>
#include <memory>

using namespace amrex;

// Adapted from amrex/Src/Base/AMReX_PlotFileDataImpl.H
// not inherited because constructor needed to be rewritten

class CheckpointFileDataImpl
{
public:
    explicit CheckpointFileDataImpl(
        std::string const& chkptfile_name,
        Vector<std::string> addl_field_names = {});

    [[nodiscard]] int spaceDim() const noexcept { return m_spacedim; }

    [[nodiscard]] Real time() const noexcept { return m_time; }

    [[nodiscard]] int finestLevel() const noexcept { return m_finest_level; }

    [[nodiscard]] const BoxArray& boxArray(int level) const noexcept
    {
        return m_ba[level];
    }

    [[nodiscard]] const DistributionMapping&
    DistributionMap(int level) const noexcept
    {
        return m_dmap[level];
    }

    void syncDistributionMap(CheckpointFileDataImpl const& src) noexcept;

    void
    syncDistributionMap(int level, CheckpointFileDataImpl const& src) noexcept;

    [[nodiscard]] Box probDomain(int level) const noexcept
    {
        return m_prob_domain[level];
    }
    [[nodiscard]] Array<Real, AMREX_SPACEDIM> probLo() const noexcept
    {
        return m_prob_lo;
    }
    [[nodiscard]] Array<Real, AMREX_SPACEDIM> probHi() const noexcept
    {
        return m_prob_hi;
    }
    [[nodiscard]] Array<Real, AMREX_SPACEDIM> cellSize(int level) const noexcept
    {
        return m_cell_size[level];
    }

    [[nodiscard]] const Vector<std::string>& varNames() const noexcept
    {
        return m_var_names;
    }

    [[nodiscard]] const Vector<std::string>& varNamesComponents() const noexcept
    {
        return m_var_names_full;
    }

    [[nodiscard]] int nComp() const noexcept { return m_ncomp; }
    [[nodiscard]] IntVect nGrowVect(int level) const noexcept
    {
        return m_ngrow[level];
    }

    MultiFab get(int level) noexcept;

private:
    std::string m_chkptfile_name;
    int m_ncomp{8};
    int m_nfields{4};
    Vector<std::string> m_var_names = {"density", "gp", "p", "velocity"};
    Vector<std::string> m_var_names_full = {
        "density", "gpx",       "gpy",       "gpz",
        "p",       "velocityx", "velocityy", "velocityz"};
    Vector<int> m_var_ncomp = {1, 3, 1, 3};
    const int m_spacedim{3};
    Real m_time;
    Real m_dt_restart;
    Real m_dt_nm1;
    Real m_dt_nm2;
    int m_finest_level, m_nlevels;
    Array<Real, AMREX_SPACEDIM> m_prob_lo{{AMREX_D_DECL(0., 0., 0.)}};
    Array<Real, AMREX_SPACEDIM> m_prob_hi{{AMREX_D_DECL(1., 1., 1.)}};
    Vector<int> m_ref_ratio;
    Vector<Box> m_prob_domain;
    int m_nstep;
    Vector<Array<Real, AMREX_SPACEDIM>> m_cell_size;
    Vector<std::string> m_mf_name;
    Vector<std::unique_ptr<VisMF>> m_vismf;
    Vector<BoxArray> m_ba;
    Vector<DistributionMapping> m_dmap;
    Vector<IntVect> m_ngrow;
};

#endif
