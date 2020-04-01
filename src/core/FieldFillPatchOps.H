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

    virtual void
    fillpatch(int lev, amrex::Real time, amrex::MultiFab& mfab) = 0;

    virtual void
    fillpatch_from_coarse(int lev, amrex::Real time, amrex::MultiFab& mfab) = 0;

    virtual void
    fillphysbc(int lev, amrex::Real time, amrex::MultiFab& mfab) = 0;
};

class FieldFillConstScalar: public FieldFillPatchOpsBase
{
public:
    FieldFillConstScalar(Field&, amrex::Real fill_val)
        : m_fill_val(fill_val)
    {}

    void fillpatch(int, amrex::Real, amrex::MultiFab& mfab) override
    {
        mfab.setVal(m_fill_val);
    }

    void fillpatch_from_coarse(
        int, amrex::Real, amrex::MultiFab& mfab) override
    {
        mfab.setVal(m_fill_val);
    }

    void fillphysbc(int, amrex::Real, amrex::MultiFab& mfab) override
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
        Field& field, amrex::AmrCore& mesh, const SimTime& time,  int probtype)
        : m_time(time), m_mesh(mesh), m_field(field), m_probtype(probtype)
    {
    }

    amrex::Vector<amrex::MultiFab*> get_mfab_vec(int lev)
    {
        const int nstates = amrex::min(m_field.num_states(), 2);
        amrex::Vector<amrex::MultiFab*> ret;
        for (int i=0; i < nstates; ++i) {
            const auto fstate = static_cast<FieldState>(i);
            ret.push_back(&m_field.state(fstate)(lev));
        }
        return ret;
    }

    void fillpatch(int lev, amrex::Real time, amrex::MultiFab& mfab) override
    {
        auto mfab_vec_lev = get_mfab_vec(lev);
        if (lev == 0) {
            amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>> physbc(
                m_mesh.Geom(lev), m_field.bcrec(),
                Functor{m_probtype, m_field.num_comp(), m_field.bc_values_device()});

            amrex::FillPatchSingleLevel(
                mfab, m_field.num_grow(), time,
                mfab_vec_lev,
                {m_time.new_time(), m_time.current_time()}, 0, 0,
                m_field.num_comp(), m_mesh.Geom(lev), physbc, 0);
        } else {
            amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>> cphysbc(
                m_mesh.Geom(lev - 1), m_field.bcrec(),
                Functor{m_probtype, m_field.num_comp(), m_field.bc_values_device()});

            amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>> fphysbc(
                m_mesh.Geom(lev), m_field.bcrec(),
                Functor{m_probtype, m_field.num_comp(), m_field.bc_values_device()});

            amrex::Interpolater* mapper = &amrex::cell_cons_interp;
            auto mfab_vec_levm1 = get_mfab_vec(lev - 1);

            amrex::FillPatchTwoLevels(
                mfab, m_field.num_grow(), time,
                mfab_vec_levm1,
                {m_time.new_time(), m_time.current_time()},
                mfab_vec_lev,
                {m_time.new_time(), m_time.current_time()}, 0, 0,
                m_field.num_comp(), m_mesh.Geom(lev - 1), m_mesh.Geom(lev),
                cphysbc, 0, fphysbc, 0, m_mesh.refRatio(lev - 1), mapper,
                m_field.bcrec(), 0);
        }
    }

    void fillpatch_from_coarse(int lev, amrex::Real time, amrex::MultiFab& mfab) override
    {
        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>> cphysbc(
            m_mesh.Geom(lev - 1), m_field.bcrec(),
            Functor{m_probtype, m_field.num_comp(), m_field.bc_values_device()});

        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>> fphysbc(
            m_mesh.Geom(lev), m_field.bcrec(),
            Functor{m_probtype, m_field.num_comp(), m_field.bc_values_device()});

        amrex::Interpolater* mapper = &amrex::cell_cons_interp;

        amrex::InterpFromCoarseLevel(
            mfab, m_field.num_grow(), time, m_field(lev), 0, 0,
            m_field.num_comp(), m_mesh.Geom(lev - 1), m_mesh.Geom(lev), cphysbc,
            0, fphysbc, 0, m_mesh.refRatio(lev - 1), mapper, m_field.bcrec(),
            0);
    }

    void fillphysbc(int lev, amrex::Real time, amrex::MultiFab& mfab) override
    {
        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>> physbc(
            m_mesh.Geom(lev), m_field.bcrec(),
            Functor{m_probtype, m_field.num_comp(), m_field.bc_values_device()});
        physbc.FillBoundary(
            mfab, 0, m_field.num_comp(), m_field.num_grow(), time, 0);
    }

protected:
    const SimTime& m_time;
    const amrex::AmrCore& m_mesh;
    Field& m_field;
    int m_probtype;
};

} // namespace amr_wind

#endif /* FIELDFILLPATCHOPS_H */
