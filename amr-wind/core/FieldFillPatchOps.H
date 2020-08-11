#ifndef FIELDFILLPATCHOPS_H
#define FIELDFILLPATCHOPS_H

#include "amr-wind/core/Field.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/core/FieldDescTypes.H"
#include "amr-wind/core/FieldUtils.H"

#include "AMReX_AmrCore.H"
#include "AMReX_MultiFab.H"
#include "AMReX_REAL.H"
#include "AMReX_PhysBCFunct.H"
#include "AMReX_FillPatchUtil.H"

/**
 *  \defgroup field_fillpatch Fillpatch operators
 *  Field fillpatch operators
 *
 *  This group documents classes and functions available for performing
 *  fillpatch operations on the entire field. This API represents a higher-level
 *  interface to the MultiFab based fillpatch operations and fill physics BC
 *  operations defined in AMReX library. The functionality available in AMR-Wind
 *  has two components: 1. the fill patch operators themselves that perform the
 *  fill operation, and 2. operator classes that define what happens during a
 *  fillpatch operation at the domain boundaries.
 *
 *  \ingroup fields
 */

namespace amr_wind {

/** Abstract interface for classes that handle fill patch and physical boundary
 *  condition operations for multi-level AMR meshes.
 *
 *  \ingroup field_fillpatch
 */
class FieldFillPatchOpsBase
{
public:
    FieldFillPatchOpsBase() = default;

    virtual ~FieldFillPatchOpsBase() = default;

    //! Implementation that handles filling patches on a single level as well as
    //! across a coarse-fine interface.
    virtual void fillpatch(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost) = 0;

    //! Implementation that handles filling patches from a coarse to fine level
    virtual void fillpatch_from_coarse(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost) = 0;

    //! Implementation that handles filling physical boundary conditions
    virtual void fillphysbc(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost) = 0;
};

/** Implementation that just fills a constant value on newly created grids
 *
 *  \ingroup field_fillpatch
 */
class FieldFillConstScalar : public FieldFillPatchOpsBase
{
public:
    FieldFillConstScalar(Field&, amrex::Real fill_val) : m_fill_val(fill_val) {}

    void fillpatch(
        int, amrex::Real, amrex::MultiFab& mfab, const amrex::IntVect&) override
    {
        mfab.setVal(m_fill_val);
    }

    void fillpatch_from_coarse(
        int, amrex::Real, amrex::MultiFab& mfab, const amrex::IntVect&) override
    {
        mfab.setVal(m_fill_val);
    }

    void fillphysbc(
        int, amrex::Real, amrex::MultiFab& mfab, const amrex::IntVect&) override
    {
        mfab.setVal(m_fill_val);
    }

private:
    amrex::Real m_fill_val;
};

/** Concrete implementation of filling patches and boundary conditions on
 *  multi-level AMR meshes.
 *
 *  \ingroup field_fillpatch
 */
template <typename Functor>
class FieldFillPatchOps : public FieldFillPatchOpsBase
{
public:
    /**
     *  @param field Field whose patches are filled by this instance
     *  @param mesh  The mesh instance to determine amrex::Geometry at a level
     *  @param time  Time instance
     *  @param probtype Unique identifier indicating problem (For BC operator)
     *  @param itype Interpolation type (default: conservative cell linear)
     */
    FieldFillPatchOps(
        Field& field,
        const amrex::AmrCore& mesh,
        const SimTime& time,
        int probtype,
        FieldInterpolator itype = FieldInterpolator::CellConsLinear)
        : m_time(time)
        , m_mesh(mesh)
        , m_field(field)
        , m_probtype(probtype)
        , m_mapper(field_impl::get_interpolation_operator(itype))
    {}

    /** Utility function to get a vector of multifabs for over states
     *
     *  Returns one or two multifabs at a given level depending on whether the
     *  field was registered with states. When two states are returned the first
     *  one is the Old state and the second one is the New state.
     */
    amrex::Vector<amrex::MultiFab*> get_mfab_vec(int lev)
    {
        const int nstates = amrex::min(m_field.num_time_states(), 2);
        amrex::Vector<amrex::MultiFab*> ret;

        // The states in the FieldInfo data are ordered from newest to oldest,
        // so swap the order
        for (int i = nstates - 1; i >= 0; --i) {
            const auto fstate = static_cast<FieldState>(i);
            ret.push_back(&m_field.state(fstate)(lev));
        }
        return ret;
    }

#if 1
    // Version that does no interpolation in time

    void fillpatch(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost) override
    {
        if (lev == 0) {
            amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>> physbc(
                m_mesh.Geom(lev), m_field.bcrec(),
                Functor{m_probtype, m_field.num_comp(),
                        m_field.bc_values_device()});

            amrex::FillPatchSingleLevel(
                mfab, nghost, time, {&m_field(lev)}, {time}, 0, 0,
                m_field.num_comp(), m_mesh.Geom(lev), physbc, 0);
        } else {
            amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>> cphysbc(
                m_mesh.Geom(lev - 1), m_field.bcrec(),
                Functor{m_probtype, m_field.num_comp(),
                        m_field.bc_values_device()});

            amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>> fphysbc(
                m_mesh.Geom(lev), m_field.bcrec(),
                Functor{m_probtype, m_field.num_comp(),
                        m_field.bc_values_device()});

            amrex::FillPatchTwoLevels(
                mfab, nghost, time, {&m_field(lev - 1)}, {time},
                {&m_field(lev)}, {time}, 0, 0, m_field.num_comp(),
                m_mesh.Geom(lev - 1), m_mesh.Geom(lev), cphysbc, 0, fphysbc, 0,
                m_mesh.refRatio(lev - 1), m_mapper, m_field.bcrec(), 0);
        }
    }

#else
    // Version that handles fields at two states and performs interpolation in
    // time.

    void fillpatch(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost) override
    {
        auto mfab_vec_lev = get_mfab_vec(lev);
        if (lev == 0) {
            amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>> physbc(
                m_mesh.Geom(lev), m_field.bcrec(),
                Functor{m_probtype, m_field.num_comp(),
                        m_field.bc_values_device()});

            amrex::FillPatchSingleLevel(
                mfab, nghost, time, mfab_vec_lev,
                {m_time.current_time(), m_time.new_time()}, 0, 0,
                m_field.num_comp(), m_mesh.Geom(lev), physbc, 0);
        } else {
            amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>> cphysbc(
                m_mesh.Geom(lev - 1), m_field.bcrec(),
                Functor{m_probtype, m_field.num_comp(),
                        m_field.bc_values_device()});

            amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>> fphysbc(
                m_mesh.Geom(lev), m_field.bcrec(),
                Functor{m_probtype, m_field.num_comp(),
                        m_field.bc_values_device()});

            auto mfab_vec_levm1 = get_mfab_vec(lev - 1);

            amrex::FillPatchTwoLevels(
                mfab, nghost, time, mfab_vec_levm1,
                {m_time.current_time(), m_time.new_time()}, mfab_vec_lev,
                {m_time.current_time(), m_time.new_time()}, 0, 0,
                m_field.num_comp(), m_mesh.Geom(lev - 1), m_mesh.Geom(lev),
                cphysbc, 0, fphysbc, 0, m_mesh.refRatio(lev - 1), m_mapper,
                m_field.bcrec(), 0);
        }
        amrex::Print() << lev << " " << m_time.current_time()
                       << " " << m_time.new_time() << std::endl;
    }
#endif

    void fillpatch_from_coarse(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost) override
    {
        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>> cphysbc(
            m_mesh.Geom(lev - 1), m_field.bcrec(),
            Functor{m_probtype, m_field.num_comp(),
                    m_field.bc_values_device()});

        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>> fphysbc(
            m_mesh.Geom(lev), m_field.bcrec(),
            Functor{m_probtype, m_field.num_comp(),
                    m_field.bc_values_device()});

        amrex::InterpFromCoarseLevel(
            mfab, nghost, time, m_field(lev - 1), 0, 0, m_field.num_comp(),
            m_mesh.Geom(lev - 1), m_mesh.Geom(lev), cphysbc, 0, fphysbc, 0,
            m_mesh.refRatio(lev - 1), m_mapper, m_field.bcrec(), 0);
    }

    void fillphysbc(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost) override
    {
        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>> physbc(
            m_mesh.Geom(lev), m_field.bcrec(),
            Functor{m_probtype, m_field.num_comp(),
                    m_field.bc_values_device()});
        physbc.FillBoundary(mfab, 0, m_field.num_comp(), nghost, time, 0);
    }

protected:
    const SimTime& m_time;
    const amrex::AmrCore& m_mesh;
    Field& m_field;
    int m_probtype;

    //! Function that handles interpolation from coarse to fine level
    amrex::Interpolater* m_mapper;
};

} // namespace amr_wind

#endif /* FIELDFILLPATCHOPS_H */
