#ifndef DIFFUSIONOPS_H
#define DIFFUSIONOPS_H

#include "amr-wind/core/MLMGOptions.H"
#include "amr-wind/equation_systems/PDETraits.H"
#include "amr-wind/equation_systems/SchemeTraits.H"
#include "amr-wind/equation_systems/PDEOps.H"
#include "amr-wind/equation_systems/PDEHelpers.H"
#include "amr-wind/diffusion/diffusion.H"

#include "AMReX_MLABecLaplacian.H"
#include "AMReX_MLTensorOp.H"
#include "AMReX_MLMG.H"
#include "AMReX_AmrCore.H"
#include "AMReX_MultiFabUtil.H"

namespace amr_wind {
namespace pde {

/** Base class for all diffusion operators
 *  \ingroup pdeop
 *
 *  This class provides the common operations for an implicit solution of a
 *  convection-diffusion equation within AMR-Wind.
 *
 *  \tparam LinOp The linear operator (see [AMREeX
 * docs](https://amrex-codes.github.io/amrex/docs_html/LinearSolvers.html))
 */
template <typename LinOp>
class DiffSolverIface
{
public:
    DiffSolverIface(
        PDEFields& fields,
        const bool has_overset,
        const std::string& prefix = "diffusion");

    virtual ~DiffSolverIface() = default;

    /** Implicit solve and update of a linear system
     *
     *  \param dt timestep size
     */
    virtual void linsys_solve(const amrex::Real dt);

    virtual void linsys_solve_impl();

protected:
    //! Sets up the linear operator (e.g., setup BCs, etc.)
    virtual void setup_operator(
        LinOp& linop,
        const amrex::Real alpha,
        const amrex::Real beta,
        const FieldState fstate);

    virtual void setup_solver(amrex::MLMG& mlmg);

    virtual void set_acoeffs(LinOp& linop, const FieldState fstate);

    template <typename L>
    void set_bcoeffs(
        L& linop,
        typename std::enable_if<
            std::is_same<L, amrex::MLTensorOp>::value>::type* = 0)
    {
        const int nlevels = m_pdefields.repo.num_active_levels();
        const auto& viscosity = m_pdefields.mueff;
        const auto& geom = m_pdefields.repo.mesh().Geom();
        const auto& mesh_fac = m_pdefields.repo.get_field("mesh_scaling_factor_cc");

        // Hold density and mesh scaling factor on faces
        std::unique_ptr<ScratchField> mesh_fac_xf, mesh_fac_yf, mesh_fac_zf;
        mesh_fac_xf = m_pdefields.repo.create_scratch_field(
            AMREX_SPACEDIM, 0, amr_wind::FieldLoc::XFACE);
        mesh_fac_yf = m_pdefields.repo.create_scratch_field(
            AMREX_SPACEDIM, 0, amr_wind::FieldLoc::YFACE);
        mesh_fac_zf = m_pdefields.repo.create_scratch_field(
            AMREX_SPACEDIM, 0, amr_wind::FieldLoc::ZFACE);

        amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>> mesh_fac_face(nlevels);

        for (int lev = 0; lev < nlevels; ++lev) {
            auto b = diffusion::average_velocity_eta_to_faces(
                geom[lev], viscosity(lev));

            // average mesh scaling factor to faces
            mesh_fac_face[lev][0] = &(*mesh_fac_xf)(lev);
            mesh_fac_face[lev][1] = &(*mesh_fac_yf)(lev);
            mesh_fac_face[lev][2] = &(*mesh_fac_zf)(lev);
            amrex::average_cellcenter_to_face(
                mesh_fac_face[lev], mesh_fac(lev), geom[lev], AMREX_SPACEDIM);

            // beta accounted for mesh mapping (x-face) = J/fac^2 * mu
            for (amrex::MFIter mfi(b[0]); mfi.isValid(); ++mfi) {
                amrex::Array4<amrex::Real> const& mu = b[0].array(mfi);
                amrex::Array4<amrex::Real const> const& fac =
                        mesh_fac_face[lev][0]->const_array(mfi);

                amrex::ParallelFor(mfi.tilebox(),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        amrex::Real det_j = fac(i,j,k,0) * fac(i,j,k,1) * fac(i,j,k,2);
                        // TODO: Store scaling factor fields on cell faces ?
                        mu(i,j,k) = mu(i,j,k) * det_j / std::pow(fac(i,j,k,0),2);
                    });
            }
            // beta accounted for mesh mapping (y-face) = J/fac^2 * mu
            for (amrex::MFIter mfi(b[1]); mfi.isValid(); ++mfi) {
                amrex::Array4<amrex::Real> const& mu = b[1].array(mfi);
                amrex::Array4<amrex::Real const> const& fac =
                        mesh_fac_face[lev][1]->const_array(mfi);

                amrex::ParallelFor(mfi.tilebox(),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        amrex::Real det_j = fac(i,j,k,0) * fac(i,j,k,1) * fac(i,j,k,2);
                        // TODO: Store scaling factor fields on cell faces ?
                        mu(i,j,k) = mu(i,j,k) * det_j / std::pow(fac(i,j,k,1),2);
                    });
            }
            // beta accounted for mesh mapping (z-face) = J/fac^2 * mu
            for (amrex::MFIter mfi(b[2]); mfi.isValid(); ++mfi) {
                amrex::Array4<amrex::Real> const& mu = b[2].array(mfi);
                amrex::Array4<amrex::Real const> const& fac =
                        mesh_fac_face[lev][2]->const_array(mfi);

                amrex::ParallelFor(mfi.tilebox(),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        amrex::Real det_j = fac(i,j,k,0) * fac(i,j,k,1) * fac(i,j,k,2);
                        // TODO: Store scaling factor fields on cell faces ?
                        mu(i,j,k) = mu(i,j,k) * det_j / std::pow(fac(i,j,k,2),2);
                    });
            }

            linop.setShearViscosity(lev, amrex::GetArrOfConstPtrs(b));
        }
    }

    template <typename L>
    void set_bcoeffs(
        L& linop,
        typename std::enable_if<
            std::is_same<L, amrex::MLABecLaplacian>::value>::type* = 0)
    {
        const int nlevels = m_pdefields.repo.num_active_levels();
        const auto& viscosity = m_pdefields.mueff;
        const auto& geom = m_pdefields.repo.mesh().Geom();
        const auto& mesh_fac = m_pdefields.repo.get_field("mesh_scaling_factor_cc");

        // Hold density and mesh scaling factor on faces
        std::unique_ptr<ScratchField> mesh_fac_xf, mesh_fac_yf, mesh_fac_zf;
        mesh_fac_xf = m_pdefields.repo.create_scratch_field(
            AMREX_SPACEDIM, 0, amr_wind::FieldLoc::XFACE);
        mesh_fac_yf = m_pdefields.repo.create_scratch_field(
            AMREX_SPACEDIM, 0, amr_wind::FieldLoc::YFACE);
        mesh_fac_zf = m_pdefields.repo.create_scratch_field(
            AMREX_SPACEDIM, 0, amr_wind::FieldLoc::ZFACE);

        amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>> mesh_fac_face(nlevels);

        for (int lev = 0; lev < nlevels; ++lev) {
            auto b = diffusion::average_velocity_eta_to_faces(
                geom[lev], viscosity(lev));

            // average mesh scaling factor to faces
            mesh_fac_face[lev][0] = &(*mesh_fac_xf)(lev);
            mesh_fac_face[lev][1] = &(*mesh_fac_yf)(lev);
            mesh_fac_face[lev][2] = &(*mesh_fac_zf)(lev);
            amrex::average_cellcenter_to_face(
                mesh_fac_face[lev], mesh_fac(lev), geom[lev], AMREX_SPACEDIM);

            // beta accounted for mesh mapping (x-face) = J/fac^2 * mu
            for (amrex::MFIter mfi(b[0]); mfi.isValid(); ++mfi) {
                amrex::Array4<amrex::Real> const& mu = b[0].array(mfi);
                amrex::Array4<amrex::Real const> const& fac =
                        mesh_fac_face[lev][0]->const_array(mfi);

                amrex::ParallelFor(mfi.tilebox(),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        amrex::Real det_j = fac(i,j,k,0) * fac(i,j,k,1) * fac(i,j,k,2);
                        // TODO: Store scaling factor fields on cell faces ?
                        mu(i,j,k) = mu(i,j,k) * det_j / std::pow(fac(i,j,k,0),2);
                    });
            }
            // beta accounted for mesh mapping (y-face) = J/fac^2 * mu
            for (amrex::MFIter mfi(b[1]); mfi.isValid(); ++mfi) {
                amrex::Array4<amrex::Real> const& mu = b[1].array(mfi);
                amrex::Array4<amrex::Real const> const& fac =
                        mesh_fac_face[lev][1]->const_array(mfi);

                amrex::ParallelFor(mfi.tilebox(),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        amrex::Real det_j = fac(i,j,k,0) * fac(i,j,k,1) * fac(i,j,k,2);
                        // TODO: Store scaling factor fields on cell faces ?
                        mu(i,j,k) = mu(i,j,k) * det_j / std::pow(fac(i,j,k,1),2);
                    });
            }
            // beta accounted for mesh mapping (z-face) = J/fac^2 * mu
            for (amrex::MFIter mfi(b[2]); mfi.isValid(); ++mfi) {
                amrex::Array4<amrex::Real> const& mu = b[2].array(mfi);
                amrex::Array4<amrex::Real const> const& fac =
                        mesh_fac_face[lev][2]->const_array(mfi);

                amrex::ParallelFor(mfi.tilebox(),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        amrex::Real det_j = fac(i,j,k,0) * fac(i,j,k,1) * fac(i,j,k,2);
                        // TODO: Store scaling factor fields on cell faces ?
                        mu(i,j,k) = mu(i,j,k) * det_j / std::pow(fac(i,j,k,2),2);
                    });
            }

            linop.setBCoeffs(lev, amrex::GetArrOfConstPtrs(b));
        }
    }

    PDEFields& m_pdefields;
    Field& m_density;

    MLMGOptions m_options;

    std::unique_ptr<LinOp> m_solver;
    std::unique_ptr<LinOp> m_applier;
};

/** Diffusion operator for scalar transport equations
 *  \ingroup pdeop
 */
template <typename PDE, typename Scheme>
struct DiffusionOp<
    PDE,
    Scheme,
    typename std::enable_if<std::is_base_of<ScalarTransport, PDE>::value>::type>
    : public DiffSolverIface<typename PDE::MLDiffOp>
{
    static_assert(
        PDE::ndim == 1, "DiffusionOp invoked for non-scalar PDE type");
    static_assert(
        std::is_same<typename PDE::MLDiffOp, amrex::MLABecLaplacian>::value,
        "Invalid linear operator for scalar diffusion operator");

    DiffusionOp(PDEFields& fields, const bool has_overset)
        : DiffSolverIface<typename PDE::MLDiffOp>(fields, has_overset)
    {
        this->m_solver->setDomainBC(
            diffusion::get_diffuse_scalar_bc(
                this->m_pdefields.field, amrex::Orientation::low),
            diffusion::get_diffuse_scalar_bc(
                this->m_pdefields.field, amrex::Orientation::high));
        this->m_applier->setDomainBC(
            diffusion::get_diffuse_scalar_bc(
                this->m_pdefields.field, amrex::Orientation::low),
            diffusion::get_diffuse_scalar_bc(
                this->m_pdefields.field, amrex::Orientation::high));
    }

    //! Computes the diffusion term that goes in the RHS
    void compute_diff_term(const FieldState fstate)
    {
        this->setup_operator(*this->m_applier, 0.0, -1.0, fstate);

        auto tau_state = std::is_same<Scheme, fvm::Godunov>::value
                             ? FieldState::New
                             : fstate;
        amrex::MLMG mlmg(*this->m_applier);
        mlmg.apply(
            this->m_pdefields.diff_term.state(tau_state).vec_ptrs(),
            this->m_pdefields.field.vec_ptrs());
    }
};

} // namespace pde
} // namespace amr_wind

#endif /* DIFFUSIONOPS_H */
