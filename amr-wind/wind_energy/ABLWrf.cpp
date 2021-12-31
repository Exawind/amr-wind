#include "amr-wind/wind_energy/ABLWrf.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/utilities/linear_interpolation.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_Print.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

namespace {

//! Return closest index (from lower) of value in vector
AMREX_FORCE_INLINE int
closest_index(const amrex::Vector<amrex::Real>& vec, const amrex::Real value)
{
    auto const it = std::upper_bound(vec.begin(), vec.end(), value);
    AMREX_ALWAYS_ASSERT(it != vec.end());

    const int idx = std::distance(vec.begin(), it);
    return std::max(idx - 1, 0);
}
} // namespace

ABLWrfForcing::ABLWrfForcing(const CFDSim& sim, const std::string identifier)
    : m_time(sim.time()), m_mesh(sim.mesh())
{
    amrex::Print() << "Constructing " << identifier << " object" << std::endl;

    amrex::ParmParse pp(identifier);
    pp.query("forcing_scheme", m_forcing_scheme);
    pp.query("control_gain", m_gain_coeff);
    amrex::Print() << "  forcing_scheme : " << m_forcing_scheme << std::endl;
    amrex::Print() << "  control_gain   : " << m_gain_coeff << std::endl;

    if(!pp.query("forcing_transition", m_forcing_transition)) {
        amrex::Print() << "  using full profile assimilation by default" << std::endl;
        m_forcing_transition = "none";
    }

    if (amrex::toLower(m_forcing_scheme) == "indirect")
    {
        if (amrex::toLower(m_forcing_transition) == "none")
        {
            if (pp.queryarr("weighting_heights", m_weighting_heights)) {
                pp.getarr("weighting_values", m_weighting_values);
                amrex::Print() << "  given weighting profile" << std::endl;
                for(int i=0; i < m_weighting_heights.size(); ++i) {
                    amrex::Print() << "  " << m_weighting_heights[i] << " " << m_weighting_values[i] << std::endl;
                }
                AMREX_ALWAYS_ASSERT(m_weighting_heights.size() == m_weighting_values.size());

            } else {
                // default is to have uniform weighting throughout
                amrex::Print() << "  setting default weighting" << std::endl;
                amrex::Real zmin = m_mesh.Geom(0).ProbLo(m_axis);
                amrex::Real zmax = m_mesh.Geom(0).ProbHi(m_axis);
                m_weighting_heights = {zmin,zmax};
                m_weighting_values = {1.0,1.0};
            }
        }
        else // weightings will be automatically set based on forcing transition
        {
            pp.get("transition_thickness", m_transition_thickness); // constant, required
            if(pp.query("constant_transition_height", m_transition_height)) {
                // set weighting profile
                setTransitionWeighting();
            } else {
                // expect to read transition_height history in netCDF input file
                // weighting profile will be updated at every step
                m_update_transition_height = true;
            }
        }

        if (pp.query("normalize_by_zmax",m_norm_zmax) && (m_norm_zmax != 0))
        {
            amrex::Real zmax = m_mesh.Geom(0).ProbHi(m_axis);
            m_scaleFact = 1.0 / zmax;
            amrex::Print() << "  set scaling factor to " << m_scaleFact << std::endl;
        }

    } // if forcing scheme is "indirect"
}

void ABLWrfForcing::setTransitionWeighting()
{
    amrex::Real zmin = m_mesh.Geom(0).ProbLo(m_axis);
    amrex::Real zmax = m_mesh.Geom(0).ProbHi(m_axis);
    m_weighting_heights = {
        zmin,
        m_transition_height,
        m_transition_height+m_transition_thickness,
        zmax
    };
    m_weighting_values = {1.0, 1.0, 0.0, 0.0};
    amrex::Print() << "setting new weighting profile" << std::endl;
    for(int i=0; i < m_weighting_heights.size(); ++i) {
        amrex::Print() << "  " << m_weighting_heights[i] << " " << m_weighting_values[i] << std::endl;
    }
}

void ABLWrfForcing::updateWeights()
{
    amrex::Print() << "Updating weights" << std::endl;
    for (int i=0; i<m_nht; ++i) {
        m_W[i] = interp::linear(m_weighting_heights, m_weighting_values, m_zht[i]);
        //amrex::Print() << "  " << m_zht[i] << " " << m_W[i] << std::endl;
    }
}

void ABLWrfForcing::indirectForcingInit()
{
    if (m_W.size() == 0)
    {
        // Will be here for:
        // - full profile assimilation
        // - partial profile assim w/ constant transition height
        // - partial profile assim w/ variable transition height (1st step only)
        amrex::Print() << "Initializing indirect forcing" << std::endl;
        m_W.resize(m_nht);
        updateWeights();
    } else if (amrex::toLower(m_forcing_transition) != "none") {
        // Will be here for:
        // - partial profile assim w/ variable transition height
        amrex::Print() << "Reinitializing indirect forcing" << std::endl;
        updateWeights();
    } else {
        amrex::Print() << "Should not be reinitializing indirect forcing!" << std::endl;
        return;
    }

    amrex::Array2D<amrex::Real, 0, 3, 0, 3> zTz;

    // Generate the matrix Z^T W Z
    for (int irow = 0; irow < 4; irow++) {
        for (int icol = 0; icol < 4; icol++) {

            zTz(irow, icol) = 0.0;

            for (int iht = 0; iht < m_nht; iht++) {
                zTz(irow, icol) =
                    zTz(irow, icol) +
                    m_W[iht] * std::pow(m_zht[iht] * m_scaleFact, (icol + irow));
            }
            //amrex::Print()<< "Z^T W Z ["<<irow<<","<<icol<<"] : " << zTz(irow,icol) << std::endl;
        }
    }
    // Invert the matrix Z^T W Z
    invertMat(zTz, m_im_zTz);
}

void ABLWrfForcing::invertMat(
    const amrex::Array2D<amrex::Real, 0, 3, 0, 3>& m,
    amrex::Array2D<amrex::Real, 0, 3, 0, 3>& im)
{
    amrex::Real A2323 = m(2, 2) * m(3, 3) - m(2, 3) * m(3, 2);
    amrex::Real A1323 = m(2, 1) * m(3, 3) - m(2, 3) * m(3, 1);
    amrex::Real A1223 = m(2, 1) * m(3, 2) - m(2, 2) * m(3, 1);
    amrex::Real A0323 = m(2, 0) * m(3, 3) - m(2, 3) * m(3, 0);
    amrex::Real A0223 = m(2, 0) * m(3, 2) - m(2, 2) * m(3, 0);
    amrex::Real A0123 = m(2, 0) * m(3, 1) - m(2, 1) * m(3, 0);
    amrex::Real A2313 = m(1, 2) * m(3, 3) - m(1, 3) * m(3, 2);
    amrex::Real A1313 = m(1, 1) * m(3, 3) - m(1, 3) * m(3, 1);
    amrex::Real A1213 = m(1, 1) * m(3, 2) - m(1, 2) * m(3, 1);
    amrex::Real A2312 = m(1, 2) * m(2, 3) - m(1, 3) * m(2, 2);
    amrex::Real A1312 = m(1, 1) * m(2, 3) - m(1, 3) * m(2, 1);
    amrex::Real A1212 = m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1);
    amrex::Real A0313 = m(1, 0) * m(3, 3) - m(1, 3) * m(3, 0);
    amrex::Real A0213 = m(1, 0) * m(3, 2) - m(1, 2) * m(3, 0);
    amrex::Real A0312 = m(1, 0) * m(2, 3) - m(1, 3) * m(2, 0);
    amrex::Real A0212 = m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0);
    amrex::Real A0113 = m(1, 0) * m(3, 1) - m(1, 1) * m(3, 0);
    amrex::Real A0112 = m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0);

    amrex::Real det =
        m(0, 0) * (m(1, 1) * A2323 - m(1, 2) * A1323 + m(1, 3) * A1223) -
        m(0, 1) * (m(1, 0) * A2323 - m(1, 2) * A0323 + m(1, 3) * A0223) +
        m(0, 2) * (m(1, 0) * A1323 - m(1, 1) * A0323 + m(1, 3) * A0123) -
        m(0, 3) * (m(1, 0) * A1223 - m(1, 1) * A0223 + m(1, 2) * A0123);
    det = 1.0 / det;

    im(0, 0) = det * (m(1, 1) * A2323 - m(1, 2) * A1323 + m(1, 3) * A1223);
    im(0, 1) = det * -(m(0, 1) * A2323 - m(0, 2) * A1323 + m(0, 3) * A1223);
    im(0, 2) = det * (m(0, 1) * A2313 - m(0, 2) * A1313 + m(0, 3) * A1213);
    im(0, 3) = det * -(m(0, 1) * A2312 - m(0, 2) * A1312 + m(0, 3) * A1212);
    im(1, 0) = det * -(m(1, 0) * A2323 - m(1, 2) * A0323 + m(1, 3) * A0223);
    im(1, 1) = det * (m(0, 0) * A2323 - m(0, 2) * A0323 + m(0, 3) * A0223);
    im(1, 2) = det * -(m(0, 0) * A2313 - m(0, 2) * A0313 + m(0, 3) * A0213);
    im(1, 3) = det * (m(0, 0) * A2312 - m(0, 2) * A0312 + m(0, 3) * A0212);
    im(2, 0) = det * (m(1, 0) * A1323 - m(1, 1) * A0323 + m(1, 3) * A0123);
    im(2, 1) = det * -(m(0, 0) * A1323 - m(0, 1) * A0323 + m(0, 3) * A0123);
    im(2, 2) = det * (m(0, 0) * A1313 - m(0, 1) * A0313 + m(0, 3) * A0113);
    im(2, 3) = det * -(m(0, 0) * A1312 - m(0, 1) * A0312 + m(0, 3) * A0112);
    im(3, 0) = det * -(m(1, 0) * A1223 - m(1, 1) * A0223 + m(1, 2) * A0123);
    im(3, 1) = det * (m(0, 0) * A1223 - m(0, 1) * A0223 + m(0, 2) * A0123);
    im(3, 2) = det * -(m(0, 0) * A1213 - m(0, 1) * A0213 + m(0, 2) * A0113);
    im(3, 3) = det * (m(0, 0) * A1212 - m(0, 1) * A0212 + m(0, 2) * A0112);
}

void ABLWrfForcing::constantForcingTransition(amrex::Vector<amrex::Real>& error) {
    // based on SOWFA6/src/ABLForcing/drivingForce/drivingForce.C
    
    // find indices
    int hLevelBlend0 = -1;
    int hLevelBlend1 = -1;
    int hLevelBlendMax = -1;
    for (int iht = 0; iht < m_nht; iht++) {
        if ((hLevelBlend1 < 0) && (m_zht[iht] >= m_transition_height)) {
            hLevelBlend1 = iht;
            hLevelBlend0 = iht - 1;
        }
        else if (m_zht[iht] >= m_transition_height + m_transition_thickness) {
            hLevelBlendMax = iht;
            break;
        }
    }

    // error checking
    if (hLevelBlend1 < 0) {
        amrex::Print() << "Note: Did not find bottom of transition layer" << std::endl;
        hLevelBlend0 = m_nht - 1;
        hLevelBlend1 = m_nht - 1;
        hLevelBlendMax = m_nht - 1;
    }
    else if (hLevelBlendMax < 0) {
        amrex::Print() << "Note: Did not find top of transition layer" << std::endl;
        hLevelBlendMax = m_nht - 1;
    }

    amrex::Print() << "Forcing transition to constant"
        << " from " << m_zht[hLevelBlend1]
        << " to " << m_zht[hLevelBlendMax] << std::endl;

    // calculate initial slope
    amrex::Real slope0 = (error[hLevelBlend1] - error[hLevelBlend0])
                       / (m_zht[hLevelBlend1] - m_zht[hLevelBlend0]);
    amrex::Real dslope = -slope0 / (m_zht[hLevelBlendMax] - m_zht[hLevelBlend1]);

    // march from hLevelBlend1 (z >= m_transition_height)
    // up to hLevelBlendMax (z >= m_transition_height + m_transition_thickness)
    // as the slope decreases linearly to 0
    for (int iht=hLevelBlend1; iht <= hLevelBlendMax; iht++)
    {
        amrex::Real slope = slope0 + dslope*(m_zht[iht] - m_zht[hLevelBlend1]);
        error[iht] = error[iht-1] + slope*(m_zht[iht] - m_zht[iht-1]);
    }

    // set the remaining levels above hLevelBlendMax to the last value
    for (int iht=hLevelBlendMax+1; iht < m_nht; iht++)
    {
        error[iht] = error[hLevelBlendMax];
    }
}

void ABLWrfForcing::blendForcings(const amrex::Vector<amrex::Real> lower, //W=1
                                  const amrex::Vector<amrex::Real> upper, //W=0
                                  amrex::Vector<amrex::Real>& error)
{
    amrex::Print() << "Blending forcings" << std::endl;
    for (int iht = 0; iht < m_nht; iht++)
    {
        error[iht] = m_W[iht] * lower[iht] + (1.0-m_W[iht]) * upper[iht];
    }
}

ABLWRFfile::ABLWRFfile(const std::string filewrf)
    : m_wrf_filename(filewrf)
{

    auto ncf = ncutils::NCFile::open_par(
        m_wrf_filename, NC_NOWRITE | NC_NETCDF4 | NC_MPIIO,
        amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL);

    m_wrf_nheight = ncf.dim("nheight").len();
    m_wrf_ntime = ncf.dim("ntime").len();
    
    m_wrf_height.resize(m_wrf_nheight);
    m_wrf_time.resize(m_wrf_ntime);

    ncf.var("heights").get(m_wrf_height.data());
    ncf.var("times").get(m_wrf_time.data());

    m_wrf_u.resize(m_wrf_nheight * m_wrf_ntime);
    m_wrf_v.resize(m_wrf_nheight * m_wrf_ntime);
    m_wrf_temp.resize(m_wrf_nheight * m_wrf_ntime);
    m_wrf_tflux.resize(m_wrf_ntime);

    ncf.var("wrf_momentum_u").get(m_wrf_u.data());
    ncf.var("wrf_momentum_v").get(m_wrf_v.data());
    ncf.var("wrf_temperature").get(m_wrf_temp.data());
    ncf.var("wrf_tflux").get(m_wrf_tflux.data());

    amrex::ParmParse pp("ABL");
    pp.query("WRF_tendency_forcing", m_abl_wrf_tendency);

}

const amrex::Vector<amrex::Real>& ABLWRFfile::wrf_heights() const
{
  return m_wrf_height;
}

const amrex::Vector<amrex::Real>& ABLWRFfile::wrf_times() const { return m_wrf_time; }

const amrex::Vector<amrex::Real>& ABLWRFfile::wrf_u() const { return m_wrf_u; }

const amrex::Vector<amrex::Real>& ABLWRFfile::wrf_v() const { return m_wrf_v; }

const amrex::Vector<amrex::Real>& ABLWRFfile::wrf_temp() const { return m_wrf_temp; }

const amrex::Vector<amrex::Real>& ABLWRFfile::wrf_tflux() const { return m_wrf_tflux; }

 bool ABLWRFfile::is_wrf_tendency_forcing() const {return m_abl_wrf_tendency; }

int ABLWRFfile::nheights() const { return m_wrf_nheight; }
int ABLWRFfile::times() const { return m_wrf_ntime; }

} // namespace amr_wind
