#ifndef SCHEMETRAITS_H
#define SCHEMETRAITS_H

#include <string>

namespace amr_wind {
namespace fvm {

struct Godunov
{
    static std::string scheme_name() { return "Godunov"; }

    static constexpr int num_states = 2;
    static constexpr int nghost_state = 3;
    static constexpr int nghost_src = 1;
    static constexpr int nghost_mac = 1;

    static constexpr int num_diff_states = 1;
    static constexpr int num_conv_states = 1;
};

struct MOL
{
    static std::string scheme_name() { return "MOL"; }

    static constexpr int num_states = 2;
    static constexpr int nghost_state = 2;
    static constexpr int nghost_src = 0;
    static constexpr int nghost_mac = 0;

    static constexpr int num_diff_states = 2;
    static constexpr int num_conv_states = 2;
};

}
}

#endif /* SCHEMETRAITS_H */
