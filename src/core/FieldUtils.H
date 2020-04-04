#ifndef FIELDUTILS_H
#define FIELDUTILS_H

#include "FieldDescTypes.H"
#include "AMReX_MultiFab.H"

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

    // Suppress warnings when compiling with CUDA
    return amrex::IndexType::TheCellType();
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

    // Suppress warnings when compiling with CUDA
    return fname;
}

inline bool is_valid_field_name(const std::string& name)
{
    // Else make sure that our guard is not in the last 8 characters
    auto found = name.rfind("__FS_");
    if (found != std::string::npos) return false;

    return true;
}

} // namespace field_impl
} // namespace amr_wind

#endif /* FIELDUTILS_H */
