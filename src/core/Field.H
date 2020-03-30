#ifndef FIELD_H
#define FIELD_H

#include <string>
#include <memory>
#include <unordered_map>

#include "AMReX_MultiFab.H"
#include "AMReX_BCRec.H"
#include "AMReX_Gpu.H"

#include "FieldDescTypes.H"

namespace amr_wind {

class Field;
class FieldRepo;

/** Information common to a field that has multiple states
 */
struct FieldInfo
{
    FieldInfo(
        const std::string& basename,
        const int ncomp,
        const int ngrow,
        const int nstates,
        const FieldLoc floc)
        : m_basename(basename)
        , m_ncomp(ncomp)
        , m_ngrow(ngrow)
        , m_nstates(nstates)
        , m_floc(floc)
        , m_bc_values(AMREX_SPACEDIM*2, amrex::Vector<amrex::Real>(ncomp))
        , m_bc_values_dview(ncomp * AMREX_SPACEDIM * 2)
        , m_bcrec(ncomp)
        , m_bcrec_d(ncomp)
    {
    }

    ~FieldInfo() = default;

    //! Copy the BC information to device data structures
    void copy_bc_to_device() noexcept;

    //! Field name without state information
    std::string m_basename;

    //! Number of components (Scalar = 1, Vector = 2, etc.)
    int m_ncomp;

    //! Ghost cells
    amrex::IntVect m_ngrow;

    //! Number of states available for this field
    int m_nstates;

    //! Cell, node, face centered field type
    FieldLoc m_floc;

    amrex::Vector<amrex::Vector<amrex::Real>> m_bc_values;
    amrex::GpuArray<const amrex::Real*, AMREX_SPACEDIM*2> m_bc_values_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_bc_values_dview;

    amrex::Vector<amrex::BCRec> m_bcrec;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_bcrec_d;

    std::unordered_map<FieldState, Field*> m_states;
};

/** Computational field
 *
 *  This is an abstraction of a collection of MultiFabs over multiple levels,
 *  and allows operations on the field as a whole. The Field cannot be
 *  constructed directly but instead must be created using field creation
 *  functions in amr_wind::FieldRepo
 */
class Field
{
public:
    friend class FieldRepo;

    Field(const Field&) = delete;
    Field& operator=(const Field&) = delete;

    //! Name of this field (including state information)
    inline const std::string& name() const { return m_name; }

    //! Number of components for this field
    inline int num_comp() const { return m_info->m_ncomp; }

    //! Ghost cells
    inline const amrex::IntVect& num_grow() const { return m_info->m_ngrow; }

    //! Number of states available for this field
    inline int num_states() const { return m_info->m_nstates; }

    //! Cell, node, or face centered field
    inline FieldLoc field_location() const { return m_info->m_floc; }

    //! State of this field instance
    inline FieldState field_state() const { return m_state; }

    //! FieldRepo instance that manages this field
    inline FieldRepo& repo() const { return m_repo; }

    //! Return reference to the host view of BC values array
    amrex::Vector<amrex::Vector<amrex::Real>>& bc_values()
    {
        return m_info->m_bc_values;
    }

    //! Return reference to host view of BCRec array
    amrex::Vector<amrex::BCRec>& bcrec() { return m_info->m_bcrec; }

    amrex::GpuArray<const amrex::Real*, AMREX_SPACEDIM*2>& bc_values_device()
    { return m_info->m_bc_values_d; }

    amrex::Gpu::DeviceVector<amrex::BCRec>& bcrec_device()
    { return m_info->m_bcrec_d; }

    //! Copy BC data from host to device
    //!
    //! The data must have been created on host using the host arrays
    //!
    inline void copy_bc_to_device() { m_info->copy_bc_to_device(); }

    //! Return field at a different time state
    Field& state(const FieldState fstate);
    const Field& state(const FieldState fstate) const;

    //! Return MultiFab instance for a given level
    amrex::MultiFab& operator()(int lev) noexcept;
    const amrex::MultiFab& operator()(int lev) const noexcept;

    //! Return a vector of MultiFab pointers for all levels
    amrex::Vector<amrex::MultiFab*> vec_ptrs() noexcept;

    amrex::Vector<const amrex::MultiFab*> vec_const_ptrs() const noexcept;

protected:
    Field(
        FieldRepo& repo,
        const std::string& name,
        const std::shared_ptr<FieldInfo>& finfo,
        FieldState state);

    //! Reference to the FieldRepository instance
    FieldRepo& m_repo;

    //! Name of the field
    std::string m_name;

    //! Common data for all field states
    std::shared_ptr<FieldInfo> m_info;

    //! State for this field
    FieldState m_state;
};

} // namespace amr_wind

#endif /* FIELD_H */
