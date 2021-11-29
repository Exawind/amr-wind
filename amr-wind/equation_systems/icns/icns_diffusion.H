#ifndef ICNS_DIFFUSION_H
#define ICNS_DIFFUSION_H

#include "amr-wind/equation_systems/PDETraits.H"
#include "amr-wind/equation_systems/PDEOps.H"
#include "amr-wind/equation_systems/PDEHelpers.H"
#include "amr-wind/equation_systems/DiffusionOps.H"
#include "amr-wind/equation_systems/icns/icns.H"
#include "amr-wind/utilities/console_io.H"

namespace amr_wind {
namespace pde {

class ICNSDiffTensorOp : public DiffSolverIface<amrex::MLTensorOp>
{
public:
    ICNSDiffTensorOp(PDEFields& fields, const bool has_overset)
        : DiffSolverIface<amrex::MLTensorOp>(fields, has_overset)
    {
        this->m_solver->setDomainBC(
            diffusion::get_diffuse_tensor_bc(
                this->m_pdefields.field, amrex::Orientation::low),
            diffusion::get_diffuse_tensor_bc(
                this->m_pdefields.field, amrex::Orientation::high));
        this->m_applier->setDomainBC(
            diffusion::get_diffuse_tensor_bc(
                this->m_pdefields.field, amrex::Orientation::low),
            diffusion::get_diffuse_tensor_bc(
                this->m_pdefields.field, amrex::Orientation::high));
    }

    template <typename Scheme>
    void compute_diff_term(const FieldState fstate)
    {
        this->setup_operator(*this->m_applier, 0.0, -1.0, fstate);

        auto tau_state = std::is_same<Scheme, fvm::Godunov>::value
                             ? FieldState::New
                             : fstate;
        auto& divtau = this->m_pdefields.diff_term.state(tau_state);

        amrex::MLMG mlmg(*this->m_applier);
        mlmg.apply(divtau.vec_ptrs(), this->m_pdefields.field.vec_ptrs());

        const auto& repo = this->m_pdefields.repo;
        const int nlevels = repo.num_active_levels();
        const auto& density = m_density.state(fstate);

        for (int lev = 0; lev < nlevels; ++lev) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(divtau(lev), amrex::TilingIfNotGPU());
                 mfi.isValid(); ++mfi) {
                const auto& bx = mfi.tilebox();
                const auto& divtau_arr = divtau(lev).array(mfi);
                const auto& rho_arr = density(lev).const_array(mfi);

                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        amrex::Real rhoinv = 1.0 / rho_arr(i, j, k);
                        divtau_arr(i, j, k, 0) *= rhoinv;
                        divtau_arr(i, j, k, 1) *= rhoinv;
                        divtau_arr(i, j, k, 2) *= rhoinv;
                    });
            }
        }
    }
};

class ICNSDiffScalarOp
{
public:
    ICNSDiffScalarOp(
        PDEFields& fields,
        const bool has_overset,
        const std::string& prefix = "diffusion")
        : m_pdefields(fields)
        , m_density(fields.repo.get_field("density"))
        , m_options(prefix, m_pdefields.field.name() + "_" + prefix)
    {
        amrex::LPInfo isolve = m_options.lpinfo();
        amrex::LPInfo iapply;

        iapply.setMaxCoarseningLevel(0);
        const auto& mesh = m_pdefields.repo.mesh();

        const auto& bclo = diffusion::get_diffuse_tensor_bc(
            m_pdefields.field, amrex::Orientation::low);
        const auto& bchi = diffusion::get_diffuse_tensor_bc(
            m_pdefields.field, amrex::Orientation::high);

        if (!has_overset) {
            m_solver_scalar.reset(new amrex::MLABecLaplacian(
                mesh.Geom(0, mesh.finestLevel()),
                mesh.boxArray(0, mesh.finestLevel()),
                mesh.DistributionMap(0, mesh.finestLevel()), isolve, {},
                AMREX_SPACEDIM));
            m_applier_scalar.reset(new amrex::MLABecLaplacian(
                mesh.Geom(0, mesh.finestLevel()),
                mesh.boxArray(0, mesh.finestLevel()),
                mesh.DistributionMap(0, mesh.finestLevel()), iapply, {},
                AMREX_SPACEDIM));
        } else {
            auto imask =
                fields.repo.get_int_field("mask_cell").vec_const_ptrs();
            m_solver_scalar.reset(new amrex::MLABecLaplacian(
                mesh.Geom(0, mesh.finestLevel()),
                mesh.boxArray(0, mesh.finestLevel()),
                mesh.DistributionMap(0, mesh.finestLevel()), imask, isolve, {},
                AMREX_SPACEDIM));
            m_applier_scalar.reset(new amrex::MLABecLaplacian(
                mesh.Geom(0, mesh.finestLevel()),
                mesh.boxArray(0, mesh.finestLevel()),
                mesh.DistributionMap(0, mesh.finestLevel()), imask, iapply, {},
                AMREX_SPACEDIM));
        }

        m_solver_scalar->setMaxOrder(m_options.max_order);
        m_applier_scalar->setMaxOrder(m_options.max_order);

        m_solver_scalar->setDomainBC(bclo, bchi);
        m_applier_scalar->setDomainBC(bclo, bchi);
    }

    template <typename Scheme>
    void compute_diff_term(const FieldState fstate)
    {
        auto tau_state = std::is_same<Scheme, fvm::Godunov>::value
                             ? FieldState::New
                             : fstate;
        const auto& repo = m_pdefields.repo;
        const int nlevels = repo.num_active_levels();
        const auto& geom = repo.mesh().Geom();

        auto& divtau = m_pdefields.diff_term.state(tau_state);
        const auto& density = m_density.state(fstate);
        const auto& viscosity = m_pdefields.mueff;
        const auto& mesh_fac = repo.get_field("mesh_scaling_factor_cc");

        const amrex::Real alpha = 0.0;
        const amrex::Real beta = -1.0;
        m_applier_scalar->setScalars(alpha, beta);

        for (int lev = 0; lev < nlevels; ++lev) {
            m_applier_scalar->setLevelBC(lev, &m_pdefields.field(lev));
        }

        //////////////////////////A coeffs/////////////////////////////////
        // TODO: Create a permanent field for det_j updated only upon regrid
        //  scratch field for determinant
        auto rho_times_det_j = repo.create_scratch_field(
            1, m_density.num_grow()[0], FieldLoc::CELL);
        for (int lev = 0; lev < nlevels; ++lev) {
            (*rho_times_det_j)(lev).setVal(0.0);
        }

        // multiply density with detJ
        for (int lev = 0; lev < nlevels; ++lev) {
            for (amrex::MFIter mfi(density(lev)); mfi.isValid(); ++mfi) {

                const auto& rho = density(lev).array(mfi);
                const auto& fac = mesh_fac(lev).array(mfi);
                const auto& prod = (*rho_times_det_j)(lev).array(mfi);

                const auto& bx = mfi.growntilebox();
                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        amrex::Real det_j =
                            fac(i, j, k, 0) * fac(i, j, k, 1) * fac(i, j, k, 2);
                        prod(i, j, k, 0) = rho(i, j, k, 0) * det_j;
                    });
            }
        }

        // set A coeffs
        for (int lev = 0; lev < nlevels; ++lev) {
            m_applier_scalar->setACoeffs(lev, (*rho_times_det_j)(lev));
        }
        ///////////////////////////////////////////////////////////////////

        //////////////////////////B coeffs/////////////////////////////////
        // Hold mesh scaling factor on faces
        std::unique_ptr<ScratchField> mesh_fac_xf, mesh_fac_yf, mesh_fac_zf;
        mesh_fac_xf = m_pdefields.repo.create_scratch_field(
            AMREX_SPACEDIM, 0, amr_wind::FieldLoc::XFACE);
        mesh_fac_yf = m_pdefields.repo.create_scratch_field(
            AMREX_SPACEDIM, 0, amr_wind::FieldLoc::YFACE);
        mesh_fac_zf = m_pdefields.repo.create_scratch_field(
            AMREX_SPACEDIM, 0, amr_wind::FieldLoc::ZFACE);

        amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>>
            mesh_fac_face(nlevels);

        // set B coeffs
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

                amrex::ParallelFor(
                    mfi.tilebox(),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        amrex::Real det_j =
                            fac(i, j, k, 0) * fac(i, j, k, 1) * fac(i, j, k, 2);
                        // TODO: Store scaling factor fields on cell faces ?
                        mu(i, j, k) =
                            mu(i, j, k) * det_j / std::pow(fac(i, j, k, 0), 2);
                    });
            }
            // beta accounted for mesh mapping (y-face) = J/fac^2 * mu
            for (amrex::MFIter mfi(b[1]); mfi.isValid(); ++mfi) {
                amrex::Array4<amrex::Real> const& mu = b[1].array(mfi);
                amrex::Array4<amrex::Real const> const& fac =
                    mesh_fac_face[lev][1]->const_array(mfi);

                amrex::ParallelFor(
                    mfi.tilebox(),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        amrex::Real det_j =
                            fac(i, j, k, 0) * fac(i, j, k, 1) * fac(i, j, k, 2);
                        // TODO: Store scaling factor fields on cell faces ?
                        mu(i, j, k) =
                            mu(i, j, k) * det_j / std::pow(fac(i, j, k, 1), 2);
                    });
            }
            // beta accounted for mesh mapping (z-face) = J/fac^2 * mu
            for (amrex::MFIter mfi(b[2]); mfi.isValid(); ++mfi) {
                amrex::Array4<amrex::Real> const& mu = b[2].array(mfi);
                amrex::Array4<amrex::Real const> const& fac =
                    mesh_fac_face[lev][2]->const_array(mfi);

                amrex::ParallelFor(
                    mfi.tilebox(),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        amrex::Real det_j =
                            fac(i, j, k, 0) * fac(i, j, k, 1) * fac(i, j, k, 2);
                        // TODO: Store scaling factor fields on cell faces ?
                        mu(i, j, k) =
                            mu(i, j, k) * det_j / std::pow(fac(i, j, k, 2), 2);
                    });
            }

            m_applier_scalar->setBCoeffs(lev, amrex::GetArrOfConstPtrs(b));
        }
        ///////////////////////////////////////////////////////////////////

        amrex::MLMG mlmg(*m_applier_scalar);
        mlmg.apply(divtau.vec_ptrs(), m_pdefields.field.vec_ptrs());

        for (int lev = 0; lev < nlevels; ++lev) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(divtau(lev), amrex::TilingIfNotGPU());
                 mfi.isValid(); ++mfi) {
                const auto& bx = mfi.tilebox();
                const auto& divtau_arr = divtau(lev).array(mfi);
                const auto& rho_arr = density(lev).const_array(mfi);

                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        amrex::Real rhoinv = 1.0 / rho_arr(i, j, k);
                        divtau_arr(i, j, k, 0) *= rhoinv;
                        divtau_arr(i, j, k, 1) *= rhoinv;
                        divtau_arr(i, j, k, 2) *= rhoinv;
                    });
            }
        }
    }

    void linsys_solve(const amrex::Real dt)
    {
        const FieldState fstate = FieldState::New;
        auto& repo = m_pdefields.repo;
        auto& field = m_pdefields.field;
        const auto& density = m_density.state(fstate);
        const int nlevels = repo.num_active_levels();
        const int ndim = field.num_comp();
        auto rhs_ptr = repo.create_scratch_field("rhs", field.num_comp(), 0);
        const auto& viscosity = m_pdefields.mueff;
        const auto& mesh_fac = repo.get_field("mesh_scaling_factor_cc");
        const auto& geom = repo.mesh().Geom();

        const amrex::Real alpha = 1.0;
        const amrex::Real beta = dt;
        m_solver_scalar->setScalars(alpha, beta);

        for (int lev = 0; lev < nlevels; ++lev) {
            m_solver_scalar->setLevelBC(lev, &m_pdefields.field(lev));
        }

        //////////////////////////A coeffs/////////////////////////////////
        // TODO: Create a permanent field for det_j updated only upon regrid
        //  scratch field for determinant
        auto rho_times_det_j = repo.create_scratch_field(
            1, m_density.num_grow()[0], FieldLoc::CELL);
        for (int lev = 0; lev < nlevels; ++lev) {
            (*rho_times_det_j)(lev).setVal(0.0);
        }

        // multiply density with detJ
        for (int lev = 0; lev < nlevels; ++lev) {
            for (amrex::MFIter mfi(density(lev)); mfi.isValid(); ++mfi) {

                const auto& rho = density(lev).array(mfi);
                const auto& fac = mesh_fac(lev).array(mfi);
                const auto& prod = (*rho_times_det_j)(lev).array(mfi);

                const auto& bx = mfi.growntilebox();
                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        amrex::Real det_j =
                            fac(i, j, k, 0) * fac(i, j, k, 1) * fac(i, j, k, 2);
                        prod(i, j, k, 0) = rho(i, j, k, 0) * det_j;
                    });
            }
        }

        // set A coeffs
        for (int lev = 0; lev < nlevels; ++lev) {
            m_solver_scalar->setACoeffs(lev, (*rho_times_det_j)(lev));
        }
        ///////////////////////////////////////////////////////////////////

        //////////////////////////B coeffs/////////////////////////////////
        // Hold mesh scaling factor on faces
        std::unique_ptr<ScratchField> mesh_fac_xf, mesh_fac_yf, mesh_fac_zf;
        mesh_fac_xf = m_pdefields.repo.create_scratch_field(
            AMREX_SPACEDIM, 0, amr_wind::FieldLoc::XFACE);
        mesh_fac_yf = m_pdefields.repo.create_scratch_field(
            AMREX_SPACEDIM, 0, amr_wind::FieldLoc::YFACE);
        mesh_fac_zf = m_pdefields.repo.create_scratch_field(
            AMREX_SPACEDIM, 0, amr_wind::FieldLoc::ZFACE);

        amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>>
            mesh_fac_face(nlevels);

        // set B coeffs
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

                amrex::ParallelFor(
                    mfi.tilebox(),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        amrex::Real det_j =
                            fac(i, j, k, 0) * fac(i, j, k, 1) * fac(i, j, k, 2);
                        // TODO: Store scaling factor fields on cell faces ?
                        mu(i, j, k) =
                            mu(i, j, k) * det_j / std::pow(fac(i, j, k, 0), 2);
                    });
            }
            // beta accounted for mesh mapping (y-face) = J/fac^2 * mu
            for (amrex::MFIter mfi(b[1]); mfi.isValid(); ++mfi) {
                amrex::Array4<amrex::Real> const& mu = b[1].array(mfi);
                amrex::Array4<amrex::Real const> const& fac =
                    mesh_fac_face[lev][1]->const_array(mfi);

                amrex::ParallelFor(
                    mfi.tilebox(),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        amrex::Real det_j =
                            fac(i, j, k, 0) * fac(i, j, k, 1) * fac(i, j, k, 2);
                        // TODO: Store scaling factor fields on cell faces ?
                        mu(i, j, k) =
                            mu(i, j, k) * det_j / std::pow(fac(i, j, k, 1), 2);
                    });
            }
            // beta accounted for mesh mapping (z-face) = J/fac^2 * mu
            for (amrex::MFIter mfi(b[2]); mfi.isValid(); ++mfi) {
                amrex::Array4<amrex::Real> const& mu = b[2].array(mfi);
                amrex::Array4<amrex::Real const> const& fac =
                    mesh_fac_face[lev][2]->const_array(mfi);

                amrex::ParallelFor(
                    mfi.tilebox(),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        amrex::Real det_j =
                            fac(i, j, k, 0) * fac(i, j, k, 1) * fac(i, j, k, 2);
                        // TODO: Store scaling factor fields on cell faces ?
                        mu(i, j, k) =
                            mu(i, j, k) * det_j / std::pow(fac(i, j, k, 2), 2);
                    });
            }

            m_solver_scalar->setBCoeffs(lev, amrex::GetArrOfConstPtrs(b));
        }
        ///////////////////////////////////////////////////////////////////

        // Always multiply with rho since there is no diffusion term for density
        for (int lev = 0; lev < nlevels; ++lev) {
            auto& rhs = (*rhs_ptr)(lev);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(rhs, amrex::TilingIfNotGPU()); mfi.isValid();
                 ++mfi) {
                const auto& bx = mfi.tilebox();
                const auto& rhs_a = rhs.array(mfi);
                const auto& fld = field(lev).const_array(mfi);
                const auto& rho = density(lev).const_array(mfi);

                amrex::ParallelFor(
                    bx, ndim,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                        rhs_a(i, j, k, n) = rho(i, j, k) * fld(i, j, k, n);
                    });
            }
        }

        amrex::MLMG mlmg(*m_solver_scalar);
        m_options(mlmg);
        mlmg.solve(
            m_pdefields.field.vec_ptrs(), rhs_ptr->vec_const_ptrs(),
            m_options.rel_tol, m_options.abs_tol);

        io::print_mlmg_info(field.name() + "_multicomponent_solve", mlmg);
    }

protected:
    PDEFields& m_pdefields;
    Field& m_density;
    MLMGOptions m_options;

    std::unique_ptr<amrex::MLABecLaplacian> m_solver_scalar;
    std::unique_ptr<amrex::MLABecLaplacian> m_applier_scalar;
};

/** Specialization of diffusion operator for ICNS
 *  \ingroup icns
 */
template <typename Scheme>
struct DiffusionOp<ICNS, Scheme>
{
    std::unique_ptr<ICNSDiffTensorOp> m_tensor_op;
    std::unique_ptr<ICNSDiffScalarOp> m_scalar_op;

    static_assert(
        ICNS::ndim == AMREX_SPACEDIM,
        "DiffusionOp invoked for scalar PDE type");

    DiffusionOp(PDEFields& fields, const bool has_overset)
    {

        bool use_tensor_op = true;
        amrex::ParmParse pp(fields.field.name() + "_diffusion");
        pp.query("use_tensor_operator", use_tensor_op);

        if (use_tensor_op) {
            m_tensor_op =
                std::make_unique<ICNSDiffTensorOp>(fields, has_overset);
        } else {
            m_scalar_op =
                std::make_unique<ICNSDiffScalarOp>(fields, has_overset);
        }
    }

    void compute_diff_term(const FieldState fstate)
    {
        if (m_tensor_op) {
            m_tensor_op->compute_diff_term<Scheme>(fstate);
        } else {
            m_scalar_op->compute_diff_term<Scheme>(fstate);
        }
    }

    void linsys_solve(const amrex::Real dt)
    {
        if (m_tensor_op)
            m_tensor_op->linsys_solve(dt);
        else
            m_scalar_op->linsys_solve(dt);
    }
};

} // namespace pde
} // namespace amr_wind

#endif /* ICNS_DIFFUSION_H */
