#ifndef ABLBOUNDARYPLANE_H
#define ABLBOUNDARYPLANE_H

#include "amr-wind/core/Field.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_Gpu.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include <AMReX_BndryRegister.H>

namespace amr_wind {

enum struct io_mode { output, input, undefined };

struct BufferData
{
    amrex::Gpu::ManagedVector<amrex::Real> data;
    amrex::Vector<size_t> start{0, 0, 0, 0};
    amrex::Vector<size_t> count{0, 0, 0, 0};
};

/** Collection of data structures and operations for reading data
 *  \ingroup we_abl
 *
 *  This class contains the inlet data structures and operations to
 *  read and interpolate inflow data.
 */
class InletData
{
    using PlaneVector = amrex::Vector<amrex::FArrayBox>;

public:
    InletData() = default;

    void resize(const int /*size*/);

    void define_plane(const amrex::Orientation /*ori*/);

    void define_level_data(
        const amrex::Orientation /*ori*/,
        const amrex::Box& /*bx*/,
        const size_t /*nc*/);

#ifdef AMR_WIND_USE_NETCDF
    void read_data(
        ncutils::NCGroup&,
        const amrex::Orientation,
        const int,
        const Field*,
        const amrex::Real,
        const amrex::Vector<amrex::Real>&);
#endif

    void read_data_native(
        const amrex::OrientationIter oit,
        amrex::BndryRegister& bndry_n,
        amrex::BndryRegister& bndry_np1,
        const int lev,
        const Field* /*fld*/,
        const amrex::Real time,
        const amrex::Vector<amrex::Real>& /*times*/);

    void interpolate(const amrex::Real /*time*/);
    bool is_populated(amrex::Orientation /*ori*/) const;
    const amrex::FArrayBox&
    interpolate_data(const amrex::Orientation ori, const int lev) const
    {
        return (*m_data_interp[ori])[lev];
    }

    int& component(const int id) { return m_components[id]; }

    int component(const int id) const { return m_components.at(id); }

    int nlevels(const amrex::Orientation ori) const
    {
        return static_cast<int>((*m_data_interp[ori]).size());
    }

    amrex::Real tn() const { return m_tn; }
    amrex::Real tnp1() const { return m_tnp1; }
    amrex::Real tinterp() const { return m_tinterp; }

private:
    amrex::Vector<std::unique_ptr<PlaneVector>> m_data_n;
    amrex::Vector<std::unique_ptr<PlaneVector>> m_data_np1;
    amrex::Vector<std::unique_ptr<PlaneVector>> m_data_interp;

    //! Time for plane at n
    amrex::Real m_tn{-1.0};

    //! Time for plane at n + 1
    amrex::Real m_tnp1{-1.0};

    //! Time for plane at interpolation
    amrex::Real m_tinterp{-1.0};

    //! Map of `{variableId : component}`
    std::unordered_map<int, int> m_components;
};

/** Interface for ABL boundary plane I/O
 *  \ingroup we_abl
 *
 *  This class performs the necessary file operations to read and
 *  write boundary planes.
 *
 *  \sa ABLFillInflow
 */
class ABLBoundaryPlane
{
    static_assert(AMREX_SPACEDIM == 3, "ABL requires 3 dimensional mesh");

public:
    explicit ABLBoundaryPlane(CFDSim& /*sim*/);

    //! Execute initialization actions after mesh has been fully generated
    void post_init_actions();

    void pre_advance_work();

    void post_advance_work();

    void initialize_data();

    void write_header();

    void write_file();

    void read_header();

    void read_file();

    void populate_data(
        const int /*lev*/,
        const amrex::Real /*time*/,
        Field& /*fld*/,
        amrex::MultiFab& /*mfab*/,
        const int dcomp = 0,
        const int orig_comp = 0) const;

#ifdef AMR_WIND_USE_NETCDF
    void impl_buffer_field(
        const amrex::Box&,
        const int,
        const int,
        const amrex::GpuArray<int, 2>&,
        const amrex::IntVect&,
        const amrex::Array4<const amrex::Real>&,
        amrex::Gpu::ManagedVector<amrex::Real>&);
#endif

    bool is_initialized() const { return m_is_initialized; }

    bool box_intersects_boundary(
        const amrex::Box& /*bx*/,
        const int /*lev*/,
        const amrex::Orientation /*ori*/) const;

private:
    const amr_wind::SimTime& m_time;
    const FieldRepo& m_repo;
    const amrex::AmrCore& m_mesh;

#ifdef AMR_WIND_USE_NETCDF
    void write_data(
        const ncutils::NCGroup& grp,
        const amrex::Orientation,
        const int,
        const Field*);
#endif

    std::string m_title{"ABL boundary planes"};

    //! Normal direction for the boundary plane
    int m_write_frequency{1};

    //! Plane names
    amrex::Vector<std::string> m_plane_names{
        {"xlo", "ylo", "zlo", "xhi", "yhi", "zhi"}};

    //! IO boundary planes
    amrex::Vector<std::string> m_planes;

    //! Start outputting after this time
    amrex::Real m_out_start_time{0.0};

#ifdef AMR_WIND_USE_NETCDF
    //! NetCDF time output counter
    size_t m_out_counter{0};
#endif

    //! File name for IO
    std::string m_filename;

    //! File name for Native time file
    std::string m_time_file;

    //! Variables for IO
    amrex::Vector<std::string> m_var_names;

    //! List of fields for IO
    amrex::Vector<Field*> m_fields;

    //! NetCDF input level time
    amrex::Vector<amrex::Real> m_in_times;
    amrex::Vector<int> m_in_timesteps;

    //! Inlet data
    InletData m_in_data;

    //! IO mode
    io_mode m_io_mode{io_mode::undefined};

    //! Flag indicating if this capability is available
    bool m_is_initialized{false};

    //! controls extents on native bndry output
    const int m_in_rad = 1;
    const int m_out_rad = 1;
    const int m_extent_rad = 0;

    //! output format for bndry output
    std::string m_out_fmt{"native"};
};

} // namespace amr_wind

#endif /* ABLBOUNDARYPLANE_H */
