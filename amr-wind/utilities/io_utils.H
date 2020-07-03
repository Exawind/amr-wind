#ifndef IO_UTILS_H
#define IO_UTILS_H

#include <string>
#include "AMReX_Vector.H"

namespace amr_wind {
namespace ioutils {

/** Generate field names for multi-component fields
 *
 *  Function generates variable names based on the number of components in the
 *  field. If the field is a scalar, it just inserts the variable name. For
 *  vectors it appends "x, y, z" to the components, and for all other fields it
 *  numbers the components start with 0.
 */
inline void add_var_names(
    amrex::Vector<std::string>& vnames,
    const std::string& fname,
    const int ncomp)
{
    const amrex::Vector<std::string> comp{"x", "y", "z"};

    switch (ncomp) {
    case 1:
        vnames.push_back(fname);
        break;

    case AMREX_SPACEDIM:
        for (int i = 0; i < AMREX_SPACEDIM; ++i)
            vnames.push_back(fname + comp[i]);
        break;

    default:
        for (int i = 0; i < ncomp; ++i)
            vnames.push_back(fname + std::to_string(i));
    }
}
}
}

#endif /* IO_UTILS_H */
