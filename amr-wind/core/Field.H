#ifndef FIELD_H
#define FIELD_H

#include <string>
#include <memory>
#include <unordered_map>
#include <algorithm>

#include "AMReX_MultiFab.H"
#include "AMReX_BCRec.H"
#include "AMReX_Gpu.H"

#include "amr-wind/incflo_enums.H"
#include "amr-wind/core/FieldDescTypes.H"
#include "amr-wind/core/ViewField.H"

namespace amr_wind {

class Field;
class FieldRepo;
class FieldFillPatchOpsBase;
class FieldBCIface;
class SimTime;

/** Information common to a field that has multiple states
 *  \ingroup fields
 */
struct FieldInfo
{
    //! Maximum number of states allowed for a field
    static constexpr int max_field_states = 5;

    FieldInfo(
        std::string basename,
        const int ncomp,
        const int ngrow,
        const int nstates,
        const FieldLoc floc);

    ~FieldInfo();

    //! Check indicating whether BCs are being set to sane values
    bool bc_initialized();

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

    ///@{
    //! Boundary condition data
    amrex::GpuArray<BC, AMREX_SPACEDIM * 2> m_bc_type;
    amrex::Vector<amrex::Vector<amrex::Real>> m_bc_values;
    amrex::GpuArray<const amrex::Real*, AMREX_SPACEDIM * 2> m_bc_values_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_bc_values_dview;

    amrex::Vector<amrex::BCRec> m_bcrec;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_bcrec_d;
    ///@}

    //! Vector holding references to available states for this field
    amrex::Vector<Field*> m_states;

    //! Function that handles filling patch and physics BC data for this field
    std::unique_ptr<FieldFillPatchOpsBase> m_fillpatch_op;

    //! Custom boundary condition actions for this field
    amrex::Vector<std::unique_ptr<FieldBCIface>> m_bc_func;

    //! Flag indicating whether BCs have been initialized and copied to device
    bool m_bc_copied_to_device{false};
};

/** Computational field
 *  \ingroup fields
 *
 *  This is an abstraction of a collection of MultiFabs over multiple levels,
 *  and allows operations on the field as a whole. The Field cannot be
 *  constructed directly but instead must be created using field creation
 *  functions in amr_wind::FieldRepo. This indirection is necessary to allow
 *  different parts of the code to hold a reference to the Field instance
 *  without having to worry about the changes in underlynig MultiFab data during
 *  regrid operations.
 *
 *  Field data is split into two groups: information that is unique to a
 *  particular field at a given time state, and data that is common to a given
 *  field with multiple timestates. The latter is held within
 *  amr_wind::FieldInfo object that is common to all time instances of a given
 *  Field (which are themselves Field instances). Different time instances of a
 *  Field can access other states through this FieldInfo instance, as well as
 *  access boundary conditions and fill patch operators applicable to all
 *  states.
 *
 *  As Field does not hold any data itself, but instead is just a thin wrapper
 *  to data held by amr_wind::FieldRepo, Field instance cannot be created
 *  directly. Instead the user must use amr_wind::FieldRepo::declare_field to
 *  create and initialize the field. This also ensures that the FieldInfo
 *  instance is managed properly during field creation. Once a field is created,
 *  the user can customize the behavior (e.g., registering fill patch operator,
 *  handling I/O on a per-field basis).
 */
class Field
{
public:
    friend class FieldRepo;

    Field(const Field&) = delete;
    Field& operator=(const Field&) = delete;
    ~Field();

    //! Name of this field (including state information)
    inline const std::string& name() const { return m_name; }

    //! Base name (without state information)
    inline const std::string& base_name() const { return m_info->m_basename; }

    //! Unique integer identifier for this field
    inline unsigned id() const { return m_id; }

    //! Number of components for this field
    inline int num_comp() const { return m_info->m_ncomp; }

    //! Ghost cells
    inline const amrex::IntVect& num_grow() const { return m_info->m_ngrow; }

    //! Number of exact time states available for this field
    inline int num_time_states() const { return m_info->m_nstates; }

    //! Number of states available for this field
    inline int num_states() const
    {
        int ns = static_cast<int>(std::count_if(
            m_info->m_states.begin(), m_info->m_states.end(),
            [](const Field* ff) { return (ff != nullptr); }));
        AMREX_ASSERT(num_time_states() <= ns);
        return ns;
    }

    //! Cell, node, or face centered field
    inline FieldLoc field_location() const { return m_info->m_floc; }

    //! State of this field instance
    inline FieldState field_state() const { return m_state; }

    //! FieldRepo instance that manages this field
    inline FieldRepo& repo() const { return m_repo; }

    inline bool& fillpatch_on_regrid() { return m_fillpatch_on_regrid; }
    inline bool fillpatch_on_regrid() const { return m_fillpatch_on_regrid; }

    //! Return true if the requested state exists for this field
    inline bool query_state(const FieldState fstate) const
    {
        const int fid = static_cast<int>(fstate);
        return (m_info->m_states[fid] != nullptr);
    }

    const amrex::GpuArray<BC, AMREX_SPACEDIM * 2>& bc_type() const
    {
        return m_info->m_bc_type;
    }

    amrex::GpuArray<BC, AMREX_SPACEDIM * 2>& bc_type()
    {
        return m_info->m_bc_type;
    }

    //! Return reference to the host view of BC values array
    amrex::Vector<amrex::Vector<amrex::Real>>& bc_values()
    {
        return m_info->m_bc_values;
    }

    //! Return reference to host view of BCRec array
    amrex::Vector<amrex::BCRec>& bcrec() const { return m_info->m_bcrec; }

    const amrex::GpuArray<const amrex::Real*, AMREX_SPACEDIM * 2>&
    bc_values_device() const
    {
        return m_info->m_bc_values_d;
    }

    const amrex::Gpu::DeviceVector<amrex::BCRec>& bcrec_device() const
    {
        return m_info->m_bcrec_d;
    }

    //! Copy BC data from host to device
    //!
    //! The data must have been created on host using the host arrays
    //!
    inline void copy_bc_to_device() { m_info->copy_bc_to_device(); }

    //! Return flag indicating whether BCs have been initialized
    inline bool bc_initialized() const { return m_info->m_bc_copied_to_device; }

    //! Return a flag indicating whether a fillpatch Op has been registered
    inline bool has_fillpatch_op() const
    {
        return static_cast<bool>(m_info->m_fillpatch_op);
    }

    /** Setup default BC conditions for fillpatch operations
     *
     *  This method initializes the necessary BC data on the field so that a
     *  fillpatch call will interpolations on periodic faces and extrapolation
     *  on all the other boundaries. This method is useful for use with derived
     *  quantities where we still require some fill patch operations.
     */
    void set_default_fillpatch_bc(
        const SimTime& time,
        const amrex::BCType::mathematicalBndryTypes bctype =
            amrex::BCType::hoextrap) noexcept;

    /** Map field to the uniform mesh
     *
     *  This method transforms the field to the uniform mesh based on
     *  the fields mesh_scaling_factor_cc or mesh_scaling_factor_nd
     *  depending on whether the field is cell-centered or node-centered,
     *  respectively.
     */
    void to_uniform_space() noexcept;

    /** Map field to the stretched mesh
     *
     *  This method transforms the field to the stretched mesh based on
     *  the fields mesh_scaling_factor_cc or mesh_scaling_factor_nd
     *  depending on whether the field is cell-centered or node-centered,
     *  respectively.
     */
    void to_stretched_space() noexcept;

    //! Return field at a different time state
    Field& state(const FieldState fstate);
    const Field& state(const FieldState fstate) const;

    //! Return MultiFab instance for a given level
    amrex::MultiFab& operator()(int lev) noexcept;
    const amrex::MultiFab& operator()(int lev) const noexcept;

    //! Return a vector of MultiFab pointers for all levels
    amrex::Vector<amrex::MultiFab*> vec_ptrs() noexcept;

    //! Return vector of `const MultiFab*` for all levels
    amrex::Vector<const amrex::MultiFab*> vec_const_ptrs() const noexcept;

    //! Advance timestep for fields with multiple states
    void advance_states() noexcept;

    //! Copy a user-specified "from_state" to "to_state"
    void copy_state(FieldState to_state, FieldState from_state) noexcept;

    //! Create a new time state after the field has been created
    Field& create_state(const FieldState fstate) noexcept;

    // NOTE: We break naming rules for methods here to provide 1-1 mapping with
    // the underlying MultiFab method names.

    //! Set the field to a constant value at all levels
    void setVal(amrex::Real value) noexcept;

    /** Set specific components of a multi-component field to a given value
     *
     *  \param value Constant value for initialization
     *  \param start_comp Starting component index (0-based)
     *  \param num_comp Number of components from starting component to
     * initialize \param num_ghost Number of ghost cells to include
     */
    void setVal(
        amrex::Real value,
        int start_comp,
        int num_comp = 1,
        int nghost = 0) noexcept;

    /** Set components to those in the vector for a multi-component field
     */
    void
    setVal(const amrex::Vector<amrex::Real>& values, int nghost = 0) noexcept;

    /** Return a sub-view of the Field instance
     */
    ViewField<Field> subview(const int scomp = 0, const int ncomp = 1)
    {
        return {*this, scomp, ncomp};
    }

    /** Register a custom fillpatch class instance for this field
     *
     *  The class must be a subclass of amr_wind::FieldFillPatchOpsBase
     */
    template <typename T, class... Args>
    void register_fill_patch_op(Args&&... args)
    {
        m_info->m_fillpatch_op.reset(new T(*this, std::forward<Args>(args)...));
    }

    /** Register a custom boundary conditions class
     */
    template <typename T, class... Args>
    void register_custom_bc(Args&&... args)
    {
        m_info->m_bc_func.emplace_back(
            new T(*this, std::forward<Args>(args)...));
    }

    void fillpatch(amrex::Real time) noexcept;
    void fillpatch(amrex::Real time, amrex::IntVect ng) noexcept;
    void fillpatch_sibling_fields(
        amrex::Real time,
        amrex::IntVect ng,
        amrex::Array<Field*, AMREX_SPACEDIM>& fields) const noexcept;

    void fillphysbc(amrex::Real time) noexcept;
    void fillphysbc(amrex::Real time, amrex::IntVect ng) noexcept;

    void apply_bc_funcs(const FieldState rho_state) noexcept;

    void fillpatch(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost) noexcept;

    void fillpatch_from_coarse(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost) noexcept;

    void fillphysbc(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost) noexcept;

    void set_inflow(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost) noexcept;

    inline void fillpatch(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const int nghost) noexcept
    {
        fillpatch(lev, time, mfab, amrex::IntVect(nghost));
    }

    inline void fillpatch_from_coarse(
        int lev, amrex::Real time, amrex::MultiFab& mfab, int nghost) noexcept
    {
        fillpatch_from_coarse(lev, time, mfab, amrex::IntVect(nghost));
    }

    inline void fillphysbc(
        int lev, amrex::Real time, amrex::MultiFab& mfab, int nghost) noexcept
    {
        fillphysbc(lev, time, mfab, amrex::IntVect(nghost));
    }

    inline void set_inflow(
        int lev, amrex::Real time, amrex::MultiFab& mfab, int nghost) noexcept
    {
        set_inflow(lev, time, mfab, amrex::IntVect(nghost));
    }

    inline bool& in_uniform_space() { return m_mesh_mapped; }
    inline bool in_uniform_space() const { return m_mesh_mapped; }

protected:
    Field(
        FieldRepo& repo,
        std::string name,
        std::shared_ptr<FieldInfo> finfo,
        const unsigned fid,
        FieldState state);

    //! Reference to the FieldRepository instance
    FieldRepo& m_repo;

    //! Name of the field
    const std::string m_name;

    //! Common data for all field states
    std::shared_ptr<FieldInfo> m_info;

    //! Unique ID for this field
    const unsigned m_id;

    //! State for this field
    const FieldState m_state;

    //! Flag indicating whether fill patch operation must be performed for this
    //! field during regrid
    bool m_fillpatch_on_regrid{false};

    //! Flag to track mesh mapping (to uniform space) of field
    bool m_mesh_mapped{false};
};

} // namespace amr_wind

#endif /* FIELD_H */
