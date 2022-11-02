#ifndef MULTIPARSER_H
#define MULTIPARSER_H

#include "amr-wind/core/vs/vector_space.H"
#include "AMReX_ParmParse.H"

namespace amr_wind::utils {

/** Utility to parse inputs from multiple namespaces.
 *
 *  This class is a thin-wrapper around AMReX ParmParse utility to query user
 *  inputs from two different namespaces (default and specialized). This allows
 *  the user to consolidate common options within a default section and then
 *  override specific values in the specialization section.
 *
 */
class MultiParser
{
public:
    MultiParser(const std::string& default_prefix, const std::string& prefix)
        : pp_default(default_prefix), pp(prefix)
    {}

    //! Return the ParmParse instance for the default namespace
    const amrex::ParmParse& default_params() const { return pp_default; }

    //! Return the ParmParse instance for the specialized namespace
    const amrex::ParmParse& params() const { return pp; }

    //! Check if the keyword is present in either namespace
    bool contains(const std::string& name) const
    {
        return pp.contains(name.c_str()) || pp_default.contains(name.c_str());
    }

    /** Get a vector from the input file
     *
     *  \param [in] name Keyword to search for
     *  \param [out] value Vector value
     */
    void get(const std::string& name, vs::Vector& value) const
    {
        amrex::Vector<vs::Vector::value_type> val;
        getarr(name, val);
        AMREX_ALWAYS_ASSERT(val.size() == AMREX_SPACEDIM);
        value.x() = val[0];
        value.y() = val[1];
        value.z() = val[2];
    }

    /** Query and return if a vector exists in the input file
     *
     *  \param [in] name Keyword to search for
     *  \param [out] value Vector value
     */
    void query(const std::string& name, vs::Vector& value) const
    {
        amrex::Vector<vs::Vector::value_type> val;
        queryarr(name, val);
        if (!val.empty()) {
            AMREX_ALWAYS_ASSERT(val.size() == AMREX_SPACEDIM);
            value.x() = val[0];
            value.y() = val[1];
            value.z() = val[2];
        }
    }

    //! Get the value for the keyword entry from either namespace
    template <typename T>
    void get(const std::string& name, T& value) const
    {
        if (pp.contains(name.c_str())) {
            pp.get(name.c_str(), value);
        } else {
            pp_default.get(name.c_str(), value);
        }
    }

    //! Get a vector of values for the given keyword entry from either namespace
    template <typename T>
    void getarr(const std::string& name, T& value) const
    {
        if (pp.contains(name.c_str())) {
            pp.getarr(name.c_str(), value);
        } else {
            pp_default.getarr(name.c_str(), value);
        }
    }

    //! Query the value for the keyword entry from either namespace
    template <typename T>
    void query(const std::string& name, T& value) const
    {
        pp_default.query(name.c_str(), value);
        pp.query(name.c_str(), value);
    }

    //! Query a vector of values for the given keyword entry from either
    //! namespace
    template <typename T>
    void queryarr(const std::string& name, T& value) const
    {
        pp_default.queryarr(name.c_str(), value);
        pp.queryarr(name.c_str(), value);
    }

private:
    amrex::ParmParse pp_default;
    amrex::ParmParse pp;
};

} // namespace amr_wind::utils

#endif /* MULTIPARSER_H */
