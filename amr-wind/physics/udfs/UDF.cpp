#include "amr-wind/physics/udfs/UDF.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/physics/udfs/LinearProfile.H"
#include "amr-wind/physics/udfs/PowerLawProfile.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace udf {

ConstValue::ConstValue(Field& fld) : m_field(fld)
{
    amrex::ParmParse pp_old("incflo");
    amrex::ParmParse pp(identifier() + "." + m_field.name());
    pp_old.queryarr(m_field.name().c_str(), m_value);
    pp.queryarr("value", m_value);

    if (fld.num_comp() != m_value.size()) {
        amrex::Abort("UDF: Invalid value for field: " + m_field.name());
    }
}

void ConstValue::operator()(int level, const amrex::Geometry&)
{
    auto& mfab = m_field(level);
    for (int i = 0; i < m_field.num_comp(); ++i) mfab.setVal(m_value[i], i, 1);
}

template <typename T>
UDFImpl<T>::UDFImpl(Field& fld) : m_field(fld), m_op(fld)
{}

template <typename T>
void UDFImpl<T>::operator()(int level, const amrex::Geometry& geom)
{
    auto& mfab = m_field(level);
    const auto& geomData = geom.data();

    const amrex::Real time = 0.0;
    const auto ncomp = m_field.num_comp();
    const auto& dop = m_op.device_instance();
    for (amrex::MFIter mfi(mfab); mfi.isValid(); ++mfi) {
        const auto& bx = mfi.tilebox();
        const auto& marr = mfab.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            for (int n = 0; n < ncomp; ++n) {
                dop({i, j, k}, marr, geomData, time, {}, n);
            }
        });
    }
}

template class UDFImpl<LinearProfile>;
template class UDFImpl<PowerLawProfile>;

} // namespace udf
} // namespace amr_wind
