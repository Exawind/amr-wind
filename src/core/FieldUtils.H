#ifndef FIELDUTILS_H
#define FIELDUTILS_H

#include "FieldDescTypes.H"
#include "AMReX_MultiFab.H"
#include "AMReX_Interpolater.H"

namespace amr_wind {
namespace field_impl {

inline amrex::IndexType index_type(const FieldLoc floc)
{
    switch (floc) {
    case FieldLoc::CELL:
        return amrex::IndexType::TheCellType();

    case FieldLoc::NODE:
        return amrex::IndexType::TheNodeType();

    case FieldLoc::XFACE:
        return amrex::IndexType(amrex::IntVect::TheDimensionVector(0));

    case FieldLoc::YFACE:
        return amrex::IndexType(amrex::IntVect::TheDimensionVector(1));

    case FieldLoc::ZFACE:
        return amrex::IndexType(amrex::IntVect::TheDimensionVector(2));
    }

#ifdef AMREX_USE_GPU
    // Suppress warnings when compiling with CUDA
    return amrex::IndexType::TheCellType();
#endif
}

inline std::string
field_name_with_state(const std::string& fname, const FieldState fstate)
{
    switch (fstate) {
    case FieldState::NP1:
        return fname;

    case FieldState::N:
        return fname + "__FS_Old";

    case FieldState::NM1:
        return fname + "__FS_NM1";

    case FieldState::NPH:
        return fname + "__FS_NPH";

    case FieldState::NMH:
        return fname + "__FS_NMH";
    }

#ifdef AMREX_USE_GPU
    // Suppress warnings when compiling with CUDA
    return fname;
#endif
}

inline bool is_valid_field_name(const std::string& name)
{
    // Else make sure that our guard is not in the last 8 characters
    auto found = name.rfind("__FS_");
    if (found != std::string::npos) return false;

    return true;
}

inline amrex::Interpolater*
get_interpolation_operator(const FieldInterpolator itype)
{
    switch (itype) {
    case FieldInterpolator::PiecewiseConstant:
        return &amrex::pc_interp;

    case FieldInterpolator::CellConsLinear:
        return &amrex::cell_cons_interp;

    case FieldInterpolator::NodeBilinear:
        return &amrex::node_bilinear_interp;
    }

    return &amrex::cell_cons_interp;
}

} // namespace field_impl
} // namespace amr_wind

#endif /* FIELDUTILS_H */
