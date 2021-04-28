#ifndef IBMODEL_H
#define IBMODEL_H

#include "amr-wind/core/Factory.H"
#include "amr-wind/immersed_boundary/IBTypes.H"
#include "amr-wind/immersed_boundary/IBOps.H"
#include "amr-wind/immersed_boundary/IBSrcOp.H"
#include "amr-wind/immersed_boundary/IBParser.H"

namespace amr_wind {

class CFDSim;

namespace ib {

/** Abstract representation of an immersed boundary source in the flowfield.
 *
 *  \ingroup immersed boundary
 *
 *  This class is intended to represent a single object (e.g., turbine, fixed
 *  wing) that exists in the computational domain. A collection of source
 *  objects (e.g., turbines in a wind farm) is modeled using the
 *  amr_wind::IB physics class that holds a collection of instances of
 *  IBModel and acts on them.
 */
class ImmersedBoundaryModel
    : public Factory<ImmersedBoundaryModel, CFDSim&, const std::string&, const int>
{
public:
    static std::string base_identifier() { return "ImmersedBoundaryModel"; }

    virtual ~ImmersedBoundaryModel() = default;

    virtual const std::string& label() const = 0;

    virtual int id() const = 0;

    virtual const IBInfo& info() const = 0;

    virtual void read_inputs(const utils::IBParser&) = 0;

    virtual void determine_influenced_procs() = 0;

    virtual void determine_root_proc(amrex::Vector<int>&) = 0;

    virtual void init_ib_source() = 0;

    virtual int num_velocity_points() const = 0;

    virtual void update_positions(VecSlice&) = 0;

    virtual void update_velocities(const VecSlice&) = 0;

    virtual void compute_forces() = 0;

    virtual void compute_source_term(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Geometry& geom) = 0;

    virtual void prepare_outputs(const std::string&) = 0;

    virtual void write_outputs() = 0;
};

/** Concrete implementation of the IBModel for different geometry types.
 *
 *  \ingroup IB
 *
 *  \tparam GeomTrait A geometry type trait that defines the behavior of the
 *  different actions (initialization, velocity sampling, force computation, and
 *  momentum source term computations) that are executed by the IBModel
 *  class through the amr_wind::ib physics driver class.
 *
 */
template <typename GeomTrait>
class IBModel : public ImmersedBoundaryModel::Register<IBModel<GeomTrait>>
{
private:
    //! Instance that holds all data related to a particular immersed boundary type
    typename GeomTrait::DataType m_data;

    //! Source computation operator
    ops::IBSrcOp<GeomTrait> m_src_op;

    //! Outputs processor
    ops::ProcessOutputsOp<GeomTrait> m_out_op;

public:
    static std::string identifier()
    {
        return GeomTrait::identifier();
    }

    /**
     *  \param sim CFD simulation environment instance
     *  \param label Unique string identifier for this instance
     *  \param id Unique integer tag for this instance
     */
    IBModel(CFDSim& sim, const std::string& label, const int id)
        : m_data(sim, label, id), m_src_op(m_data), m_out_op(m_data)
    {}

    //! Return the unique identifier (name) for this instance
    const std::string& label() const override { return m_data.info().label; }

    //! Return the unique tag (integer ID) for this instance
    int id() const override { return m_data.info().id; }

    //! Return the core info object for this actuator instance
    const IBInfo& info() const override { return m_data.info(); }

    void read_inputs(const utils::IBParser& pp) override
    {
        ops::ReadInputsOp<GeomTrait>()(m_data, pp);
        m_out_op.read_io_options(pp);
    }

    void determine_influenced_procs() override
    {
        ops::determine_influenced_procs<GeomTrait>(m_data);
    }

    void determine_root_proc(amrex::Vector<int>& ib_proc_count) override;

    int num_velocity_points() const override;

    void update_positions(VecSlice& pos) override
    {
        ops::UpdatePosOp<GeomTrait>()(m_data);
        const auto& vpos = m_data.grid().vel_pos;
        std::copy(vpos.begin(), vpos.end(), pos.begin());
    }

    void update_velocities(const VecSlice& vel) override
    {
        std::copy(vel.begin(), vel.end(), m_data.grid().vel.begin());
        ops::UpdateVelOp<GeomTrait>()(m_data);
    }

    void compute_forces() override
    {
        ops::ComputeForceOp<GeomTrait>()(m_data);
        m_src_op.setup_op();
    }

    void compute_source_term(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Geometry& geom) override
    {
        m_src_op(lev, mfi, geom);
    }

    void prepare_outputs(const std::string& out_dir) override
    {
        m_out_op.prepare_outputs(out_dir);
    }

    void write_outputs() override { m_out_op.write_outputs(); }

    void init_ib_source() override
    {
        ops::InitDataOp<GeomTrait>()(m_data);
        m_src_op.initialize();
    }
};

template <typename GeomTrait>
void IBModel<GeomTrait>::determine_root_proc(
    amrex::Vector<int>& ib_proc_count)
{
    ops::determine_root_proc<GeomTrait>(m_data, ib_proc_count);
    {
        // Sanity checks
        const auto& info = m_data.info();
        const auto& plist = info.procs;
        AMREX_ALWAYS_ASSERT(info.root_proc > -1);
        AMREX_ALWAYS_ASSERT(plist.find(info.root_proc) != plist.end());
    }
}

template <typename GeomTrait>
int IBModel<GeomTrait>::num_velocity_points() const
{
    auto& info = m_data.info();
    auto& grid = m_data.grid();

    return (info.sample_vel_in_proc) ? grid.vel.size() : 0;
}

} // namespace ib
} // namespace amr_wind

#endif /* IBMODEL_H */
