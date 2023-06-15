#ifndef SAMPLINGUTILS_H
#define SAMPLINGUTILS_H

#include <string>
#include "AMReX_Vector.H"
#include "amr-wind/core/vs/vector_space.H"
#include "amr-wind/utilities/tensor_ops.H"

namespace amr_wind::sampling::sampling_utils {

/** Utility functions available to advanced sampling techniques
 *
 *  \param key Name of the PDE
 *  \param terms List of identifiers for the source terms to be added
 */
// void inject_turbulence_src_terms(
//     const std::string& key, const amrex::Vector<std::string>& terms);

enum class NormalRule { SIGMA1, SIGMA2, SIGMA3, HALFPOWER };

vs::Vector
rotate_euler_vector(vs::Vector& axis, double& angle, vs::Vector& vec);

vs::Vector canon_rotator(const vs::Vector& angles, const vs::Vector& data);
vs::Tensor unit_projection_matrix(const vs::Vector& a);
vs::Tensor rotation_matrix(vs::Vector dst, vs::Vector src);
vs::Tensor skew_cross(vs::Vector a, vs::Vector b);
vs::Tensor scale(vs::Tensor v, double a);

void spherical_cap_quadrature(
    double gammav,
    int ntheta,
    std::vector<double> abscissae1D,
    std::vector<double> weights1D,
    std::vector<vs::Vector>& rays,
    std::vector<double>& weights);

void spherical_cap_truncated_normal(
    double gammav,
    int ntheta,
    NormalRule rule,
    std::vector<vs::Vector>& rays,
    std::vector<double>& weights);

std::pair<std::vector<double>, std::vector<double>>
truncated_normal_rule(NormalRule rule);

} // namespace amr_wind::sampling::sampling_utils

#endif /* SAMPLINGUTILS_H */
