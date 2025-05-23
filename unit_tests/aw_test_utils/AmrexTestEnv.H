#ifndef AMREXTESTENV_H
#define AMREXTESTENV_H

#include "gtest/gtest.h"
#include "AMReX.H"
#include "AMReX_ParmParse.H"

namespace amr_wind_tests {

/** Global setup/teardown test environment configuration for AMR-Wind unit tests
 *
 *  This class is registered with GoogleTest infrastructure to perform global
 *  setup/teardown tasks. The base implementation calls the amrex::Initialize
 *  and amrex::Finalize calls.
 *
 *  During the AmrexTestEnv::SetUp call, it also finalizes the amrex::ParmParse
 *  global instance so that each test can utilize a clean "input file". The user
 *  can disable this feature by passing `utest.keep_parameters=1` at the command
 *  line.
 *
 */
class AmrexTestEnv : public ::testing::Environment
{
public:
    AmrexTestEnv(int& argc, char**& argv) : m_argc(argc), m_argv(argv) {}

    ~AmrexTestEnv() override = default;

    void SetUp() override
    {
        using namespace amrex::mpidatatypes;
        amrex::Initialize(m_argc, m_argv, true, MPI_COMM_WORLD, []() {
            amrex::ParmParse pp("amrex");
            if (!(pp.contains("v") || pp.contains("verbose"))) {
                pp.add("verbose", -1);
                pp.add("v", -1);
            }

            pp.add("throw_exception", 1);
            pp.add("signal_handling", 0);
        });

        // Save managed memory flag for future use
        {
            amrex::ParmParse pp("amrex");
            pp.query("the_arena_is_managed", m_has_managed_memory);
        }

        // Call ParmParse::Finalize immediately to allow unit tests to start
        // with a clean "input file". However, allow user to override this
        // behavior through command line arguments.
        {
            amrex::ParmParse pp("utest");
            bool keep_parameters = false;
            pp.query("keep_parameters", keep_parameters);

            if (!keep_parameters) {
                amrex::ParmParse::Finalize();
            }
        }
    }

    void TearDown() override { amrex::Finalize(); }

    bool has_managed_memory() const { return m_has_managed_memory; }

protected:
    int& m_argc;
    char**& m_argv;

    bool m_has_managed_memory{true};
};

} // namespace amr_wind_tests

#endif /* AMREXTESTENV_H */
