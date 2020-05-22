#include "derive_K.H"
#include "FieldRepo.H"
namespace amr_wind {

template<typename FType>
void normalize_field(FType& Field)
{
    const auto& repo = Field.repo();
    const auto ncomp = Field.num_comp();
    
    const int nlevels = repo.num_active_levels();
    for (int lev=0; lev < nlevels; ++lev) {
        for (amrex::MFIter mfi(Field(lev)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.growntilebox(Field.num_grow());
            const auto& field_arr = Field(lev).array(mfi);

            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real mag=0;
                //Compute magnitude
                for (int icomp=0;icomp<ncomp;++icomp){
                    mag=mag+field_arr(i,j,k,icomp)*field_arr(i,j,k,icomp);
                }
                //Normalize field
                for (int icomp=0;icomp<ncomp;++icomp){
                    field_arr(i,j,k)=field_arr(i,j,k,icomp)/std::sqrt(mag);
                }

                });
        }
    }
}

template void normalize_field<Field>(Field&);
template void normalize_field<ScratchField>(ScratchField&);

}
