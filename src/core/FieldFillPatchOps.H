#ifndef FIELDFILLPATCHOPS_H
#define FIELDFILLPATCHOPS_H

#include "Field.H"
#include "SimTime.H"

#include "AMReX_AmrCore.H"
#include "AMReX_MultiFab.H"
#include "AMReX_REAL.H"
#include "AMReX_PhysBCFunct.H"
#include "AMReX_FillPatchUtil.H"

namespace amr_wind {

class FieldFillPatchOpsBase
{
public:
    FieldFillPatchOpsBase() = default;

    virtual ~FieldFillPatchOpsBase() = default;

    virtual void fillpatch(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost) = 0;

    virtual void fillpatch_from_coarse(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost) = 0;

    virtual void fillphysbc(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost) = 0;
};

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

template <typename Functor>
class FieldFillPatchOps : public FieldFillPatchOpsBase
{
public:
    FieldFillPatchOps(
        Field& field, amrex::AmrCore& mesh, const SimTime& time, int probtype)
        : m_time(time), m_mesh(mesh), m_field(field), m_probtype(probtype)
    {}

    /** Utility function to get a vector of multifabs for over states
     *
     *  Returns one or two multifabs at a given level depending on whether the
     *  field was registered with states. When two states are returned the first
     *  one is the Old state and the second one is the New state.
     */
    amrex::Vector<amrex::MultiFab*> get_mfab_vec(int lev)
    {
        const int nstates = amrex::min(m_field.num_states(), 2);
        amrex::Vector<amrex::MultiFab*> ret;

        // The states in the FieldInfo data are ordered from newest to oldest,
        // so swap the order
        for (int i = nstates - 1; i >= 0; --i) {
            const auto fstate = static_cast<FieldState>(i);
            ret.push_back(&m_field.state(fstate)(lev));
        }
        return ret;
    }

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

            amrex::Interpolater* mapper = &amrex::cell_cons_interp;
            auto mfab_vec_levm1 = get_mfab_vec(lev - 1);

            amrex::FillPatchTwoLevels(
                mfab, nghost, time, mfab_vec_levm1,
                {m_time.current_time(), m_time.new_time()}, mfab_vec_lev,
                {m_time.current_time(), m_time.new_time()}, 0, 0,
                m_field.num_comp(), m_mesh.Geom(lev - 1), m_mesh.Geom(lev),
                cphysbc, 0, fphysbc, 0, m_mesh.refRatio(lev - 1), mapper,
                m_field.bcrec(), 0);
        }
    }

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

        amrex::Interpolater* mapper = &amrex::cell_cons_interp;

        amrex::InterpFromCoarseLevel(
            mfab, nghost, time, m_field(lev - 1), 0, 0, m_field.num_comp(),
            m_mesh.Geom(lev - 1), m_mesh.Geom(lev), cphysbc, 0, fphysbc, 0,
            m_mesh.refRatio(lev - 1), mapper, m_field.bcrec(), 0);
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
};

} // namespace amr_wind

#endif /* FIELDFILLPATCHOPS_H */
