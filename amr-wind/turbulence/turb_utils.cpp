#include "turb_utils.H"

#include <set>

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace turbulence {
namespace turb_utils {

void inject_turbulence_src_terms(
    const std::string& key, const amrex::Vector<std::string>& terms)
{
    amrex::ParmParse pp(key);
    amrex::Vector<std::string> srcin;
    pp.queryarr("sources", srcin);

    // If there are no sources (usually the case) then add sources and return
    if (srcin.size() < 1) {
        pp.addarr("sources", terms);
        return;
    }

    // There are source terms defined by user... so need to combine terms
    std::set<std::string> sterm_set;
    for (const auto& sname: srcin)
        sterm_set.insert(sname);
    for (const auto& sname: terms)
        sterm_set.insert(sname);

    // Now convert the combined terms back into a vector for insertion in input dictionary
    srcin.clear();
    for (const auto& sname: sterm_set)
        srcin.emplace_back(sname);
    pp.addarr("sources", srcin);
}

}
}
}
