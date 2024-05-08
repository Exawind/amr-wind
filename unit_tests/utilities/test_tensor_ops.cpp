#include "aw_test_utils/AmrexTest.H"
#include "amr-wind/utilities/tensor_ops.H"

namespace amr_wind_tests {
namespace {
constexpr double tol = 1.0e-12;

struct TestVector
{

    TestVector()
    {
        AMREX_ALWAYS_ASSERT(m_hvec.size() == 3);
        m_dvec.resize(m_hvec.size());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_hvec.begin(), m_hvec.end(),
            m_dvec.begin());
    }

    ~TestVector() = default;

    amrex::Real* data() { return m_dvec.data(); }
    const amrex::Real* data() const { return m_dvec.data(); }
    size_t size() const { return m_hvec.size(); }

private:
    const amrex::Vector<amrex::Real> m_hvec = {1, 2, 3};
    amrex::Gpu::DeviceVector<amrex::Real> m_dvec;
};
} // namespace

TEST(TensorOps, vec_norm)
{
    const amrex::Real expected_value = 14;

    const TestVector tv;
    const amrex::Real* pvec = tv.data();
    amrex::Gpu::DeviceScalar<amrex::Real> ds(0.0);
    auto* ddata = ds.dataPtr();
    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int /*unused*/) {
        ddata[0] = amr_wind::utils::vec_norm(pvec);
    });
    EXPECT_NEAR(ds.dataValue(), expected_value, tol);
}

TEST(TensorOps, vec_mag)
{
    const amrex::Real expected_value = std::sqrt(14);

    const TestVector tv;
    const amrex::Real* pvec = tv.data();
    amrex::Gpu::DeviceScalar<amrex::Real> ds(0.0);
    auto* ddata = ds.dataPtr();
    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int /*unused*/) {
        ddata[0] = amr_wind::utils::vec_mag(pvec);
    });
    EXPECT_NEAR(ds.dataValue(), expected_value, tol);
}

TEST(TensorOps, vec_normalize)
{
    const amrex::Real mag = std::sqrt(14);
    const amrex::Vector<amrex::Real> expected_values = {
        1.0 / mag, 2.0 / mag, 3.0 / mag};

    TestVector tv;
    amrex::Real* pvec = tv.data();
    amrex::Gpu::DeviceVector<amrex::Real> ds(tv.size(), 0.0);
    auto* ddata = ds.dataPtr();
    const auto np = tv.size();
    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int /*unused*/) {
        amr_wind::utils::vec_normalize(pvec);
        for (int i = 0; i < np; i++) {
            ddata[i] = pvec[i];
        }
    });

    amrex::Vector<amrex::Real> hs(ds.size(), 0.0);
    amrex::Gpu::copy(
        amrex::Gpu::deviceToHost, ds.begin(), ds.end(), hs.begin());
    for (int i = 0; i < hs.size(); i++) {
        EXPECT_NEAR(hs[i], expected_values[i], tol);
    }
}
TEST(TensorOps, dot_prod)
{
    const amrex::Real expected_value = 14;

    const TestVector tv;
    const amrex::Real* pvec = tv.data();
    amrex::Gpu::DeviceScalar<amrex::Real> ds(0.0);
    auto* ddata = ds.dataPtr();
    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int /*unused*/) {
        ddata[0] = amr_wind::utils::dot_prod(pvec, pvec);
    });
    EXPECT_NEAR(ds.dataValue(), expected_value, tol);
}

TEST(TensorOps, cross_prod)
{
    const amrex::Vector<amrex::Real> expected_values = {0.0, 0.0, 0.0};

    const TestVector tv;
    const amrex::Real* pvec = tv.data();
    amrex::Gpu::DeviceVector<amrex::Real> ds(tv.size(), -1.0);
    auto* ddata = ds.dataPtr();
    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int /*unused*/) {
        amr_wind::utils::cross_prod(pvec, pvec, ddata);
    });

    amrex::Vector<amrex::Real> hs(ds.size(), -1.0);
    amrex::Gpu::copy(
        amrex::Gpu::deviceToHost, ds.begin(), ds.end(), hs.begin());
    for (int i = 0; i < hs.size(); i++) {
        EXPECT_NEAR(hs[i], expected_values[i], tol);
    }
}

TEST(TensorOps, transform_vec)
{
    const amrex::Vector<amrex::Real> expected_values = {14.0, -20.0, 32.0};

    const amrex::Real tmat[AMREX_SPACEDIM][AMREX_SPACEDIM] = {
        {1.0, 2.0, 3.0}, {-2.0, -3.0, -4.0}, {4.0, 5.0, 6.0}};

    const TestVector tv;
    const amrex::Real* pvec = tv.data();
    amrex::Gpu::DeviceVector<amrex::Real> ds(tv.size(), -1.0);
    auto* ddata = ds.dataPtr();
    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int /*unused*/) {
        amr_wind::utils::transform_vec(tmat, pvec, ddata);
    });

    amrex::Vector<amrex::Real> hs(ds.size(), -1.0);
    amrex::Gpu::copy(
        amrex::Gpu::deviceToHost, ds.begin(), ds.end(), hs.begin());
    for (int i = 0; i < hs.size(); i++) {
        EXPECT_NEAR(hs[i], expected_values[i], tol);
    }
}

TEST(TensorOps, inv_transform_vec)
{
    const amrex::Vector<amrex::Real> expected_values = {9.0, 11.0, 13.0};

    const amrex::Real tmat[AMREX_SPACEDIM][AMREX_SPACEDIM] = {
        {1.0, 2.0, 3.0}, {-2.0, -3.0, -4.0}, {4.0, 5.0, 6.0}};

    const TestVector tv;
    const amrex::Real* pvec = tv.data();
    amrex::Gpu::DeviceVector<amrex::Real> ds(tv.size(), -1.0);
    auto* ddata = ds.dataPtr();
    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int /*unused*/) {
        amr_wind::utils::inv_transform_vec(tmat, pvec, ddata);
    });

    amrex::Vector<amrex::Real> hs(ds.size(), -1.0);
    amrex::Gpu::copy(
        amrex::Gpu::deviceToHost, ds.begin(), ds.end(), hs.begin());
    for (int i = 0; i < hs.size(); i++) {
        EXPECT_NEAR(hs[i], expected_values[i], tol);
    }
}

} // namespace amr_wind_tests
