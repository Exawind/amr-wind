#ifndef ACTPARSER_H
#define ACTPARSER_H

#include "amr-wind/core/vs/vector_space.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace actuator {
namespace utils {

/** Utility to parse actuator inputs from multiple namespaces
 *
 *  The actuator inputs can be provided within a default namespace (e.g.,
 *  TurbineLine, TurbineDisk) or a namespace specific to the turbine (e.g.,
 *  Turbine01). This reduces duplication in input file and allows overriding
 *  values for specific turbines.
 */
class ActParser
{
public:
    ActParser(const std::string& default_prefix, const std::string& prefix)
        : pp_default(default_prefix), pp(prefix)
    {}

    const amrex::ParmParse& default_params() const { return pp_default; }
    const amrex::ParmParse& params() const { return pp; }

    //! Provide special accessor for obtaining vectors from ParmParse namespaces
    void get(const std::string& name, vs::Vector& value) const
    {
        amrex::Vector<vs::Vector::value_type> val;
        getarr(name, val);
        AMREX_ALWAYS_ASSERT(val.size() == AMREX_SPACEDIM);
        value.x() = val[0];
        value.y() = val[1];
        value.z() = val[2];
    }

    void query(const std::string& name, vs::Vector& value) const
    {
        amrex::Vector<vs::Vector::value_type> val;
        queryarr(name, val);
        if (val.size() > 0) {
            AMREX_ALWAYS_ASSERT(val.size() == AMREX_SPACEDIM);
            value.x() = val[0];
            value.y() = val[1];
            value.z() = val[2];
        }
    }

    template<typename T>
    void get(const std::string& name, T& value) const
    {
        if (pp.contains(name.c_str())) {
            pp.get(name.c_str(), value);
        } else {
            pp_default.get(name.c_str(), value);
        }
    }

    template<typename T>
    void getarr(const std::string& name, T& value) const
    {
        if (pp.contains(name.c_str())) {
            pp.getarr(name.c_str(), value);
        } else {
            pp_default.getarr(name.c_str(), value);
        }
    }

    template<typename T>
    void query(const std::string& name, T& value) const
    {
        pp_default.query(name.c_str(), value);
        pp.query(name.c_str(), value);
    }

    template<typename T>
    void queryarr(const std::string& name, T& value) const
    {
        pp_default.queryarr(name.c_str(), value);
        pp.queryarr(name.c_str(), value);
    }

private:
    amrex::ParmParse pp_default;
    amrex::ParmParse pp;
};

} // namespace utils
} // namespace actuator
} // namespace amr_wind

#endif /* ACTPARSER_H */
