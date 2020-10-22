#include "amr-wind/utilities/DerivedQuantity.H"
#include "amr-wind/core/field_ops.H"
#include "amr-wind/fvm/fvm.H"
#include "amr-wind/utilities/io_utils.H"

namespace amr_wind {
namespace derived {

struct VorticityMag : public DerivedQty::Register<VorticityMag>
{
    static constexpr int ncomp = 1;
    static std::string identifier() { return "mag_vorticity"; }

    VorticityMag(const FieldRepo& repo, const std::vector<std::string>& args)
        : m_vel(repo.get_field("velocity"))
    {
        AMREX_ALWAYS_ASSERT(args.size() == 0u);
    }

    std::string name() const override { return identifier(); }

    int num_comp() const override { return ncomp; }

    void operator()(ScratchField& fld, const int scomp = 0) override
    {
        AMREX_ASSERT(fld.num_comp() > (scomp));
        auto vort_mag = fld.subview(scomp, 1);
        fvm::vorticity_mag(vort_mag, m_vel);
    }

private:
    const Field& m_vel;
};

struct StrainRateMag : public DerivedQty::Register<StrainRateMag>
{
    static constexpr int ncomp = 1;
    static std::string identifier() { return "mag_strainrate"; }

    StrainRateMag(const FieldRepo& repo, const std::vector<std::string>& args)
        : m_vel(repo.get_field("velocity"))
    {
        AMREX_ALWAYS_ASSERT(args.size() == 0u);
    }

    std::string name() const override { return identifier(); }

    int num_comp() const override { return ncomp; }

    void operator()(ScratchField& fld, const int scomp = 0) override
    {
        AMREX_ASSERT(fld.num_comp() > (scomp));
        auto srate = fld.subview(scomp, 1);
        fvm::strainrate(srate, m_vel);
    }

private:
    const Field& m_vel;
};

struct Gradient : public DerivedQty::Register<Gradient>
{
    static std::string identifier() { return "grad"; }

    Gradient(const FieldRepo& repo, const std::vector<std::string>& args)
    {
        AMREX_ALWAYS_ASSERT(args.size() == 1u);
        m_phi = &repo.get_field(args[0]);
    }

    std::string name() const override { return "grad_" + m_phi->name(); }

    int num_comp() const override
    { return AMREX_SPACEDIM * m_phi->num_comp(); }

    void operator()(ScratchField& fld, const int scomp = 0) override
    {
        AMREX_ASSERT(fld.num_comp() >= (scomp + num_comp()));
        auto gradphi = fld.subview(scomp, num_comp());
        fvm::gradient(gradphi, *m_phi);
    }

private:
    const Field* m_phi;
};

struct Divergence : public DerivedQty::Register<Divergence>
{
    static std::string identifier() { return "div"; }

    Divergence(const FieldRepo& repo, const std::vector<std::string>& args)
    {
        AMREX_ALWAYS_ASSERT(args.size() == 1u);
        m_phi = &repo.get_field(args[0]);
        AMREX_ALWAYS_ASSERT(m_phi->num_comp() == AMREX_SPACEDIM);
    }

    std::string name() const override { return "div_" + m_phi->name(); }

    int num_comp() const override { return 1; }

    void operator()(ScratchField& fld, const int scomp = 0) override
    {
        AMREX_ASSERT(fld.num_comp() >= (scomp + num_comp()));
        auto divphi = fld.subview(scomp, num_comp());
        fvm::divergence(divphi, *m_phi);
    }

private:
    const Field* m_phi;
};

struct Laplacian : public DerivedQty::Register<Laplacian>
{
    static constexpr int ncomp = 1;
    static std::string identifier() { return "laplacian"; }

    Laplacian(const FieldRepo& repo, const std::vector<std::string>& args)
    {
        AMREX_ALWAYS_ASSERT(args.size() == 1u);
        m_phi = &repo.get_field(args[0]);
        AMREX_ALWAYS_ASSERT(m_phi->num_comp() == AMREX_SPACEDIM);
    }

    std::string name() const override { return "laplacian_" + m_phi->name(); }

    int num_comp() const override { return ncomp; }

    void operator()(ScratchField& fld, const int scomp = 0) override
    {
        AMREX_ASSERT(fld.num_comp() >= (scomp + num_comp()));
        auto lapphi = fld.subview(scomp, num_comp());
        fvm::laplacian(lapphi, *m_phi);
    }

private:
    const Field* m_phi;
};

struct FieldComponents : public DerivedQty::Register<FieldComponents>
{
    static std::string identifier() { return "components"; }

    FieldComponents(const FieldRepo& repo, const std::vector<std::string>& args)
    {
        const size_t nargs = args.size();
        AMREX_ALWAYS_ASSERT(nargs > 1u);
        m_fld = &repo.get_field(args[0]);
        AMREX_ALWAYS_ASSERT(static_cast<int>(nargs - 1) < m_fld->num_comp());

        m_ncomp = nargs - 1;
        m_comp.resize(nargs - 1);
        for (size_t i=1; i < nargs; ++i) {
            m_comp[i - 1] = std::stoi(args[i]);
            AMREX_ALWAYS_ASSERT((m_comp[i-1] >= 0) && (m_comp[i - 1] < m_fld->num_comp()));
        }
    }

    std::string name() const override { return m_fld->name(); }

    int num_comp() const override  { return m_ncomp; }

    void var_names(amrex::Vector<std::string>& plt_var_names) override
    {
        amrex::Vector<std::string> names;
        ioutils::add_var_names(names, m_fld->name(), m_fld->num_comp());
        for (auto ic: m_comp) {
            plt_var_names.push_back(names[ic]);
        }
    }

    void operator()(ScratchField& fld, const int scomp = 0) override
    {
        AMREX_ASSERT(fld.num_comp() >= (scomp + num_comp()));
        int dst_comp = scomp;
        for (auto icomp: m_comp) {
            field_ops::copy(fld, *m_fld, icomp, dst_comp, 1, 0);
            ++dst_comp;
        }
    }

private:
    const Field* m_fld;
    amrex::Vector<int> m_comp;
    int m_ncomp{0};
};

} // namespace derived
} // namespace amr_wind
