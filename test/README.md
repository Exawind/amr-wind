# incflo regression tests

Regression tests for incflo are managed using [the regression testing tools that are distributed as part of AMReX](https://github.com/AMReX-Codes/regression_testing).  
The tools are a custom set of python scripts that check the
results produced by the code in the remote repository (i.e. not in your local, working copy).  
For information on how to set up a remote repository, see the Setting Up the Tests section below.

The regression test scripts
1. "git pull" AMReX and incflo
2. For each defined test,
   * build an executable according to compile-time parameters/defines
   * run in serial or parallel
   * compare the results (typically in the form of a plotfile) with a
      "benchmark" reference solution.
3. Assemble the results of all the tests, formatted in html.
   Optionally, if any of the tests fail, an email is generated and sent
   to a specified list of recipients.

A key design feature of the regression suite is that the reference
solutions can be updated manually at any time.  This is necessary
when, for example, a bug is discovered or an algorithm change results
in improved solutions or modified error metrics.  For this reason, the
initial set of benchmarks needs to be created manually before the first
test is executed (more info on how to do this in the Setting Up the
Tests section).


## Setting Up the Tests

1. Choose/create an area/directory to run the regression tests.  It is
recommended to make a special "scratch" area for the exclusive use of
the regression tester.  This is because the testing suite may optionally
checkout specific branches or SHA1 commits of the needed repositories,
and while it always restores the repository to it's original state, it's
best to let the tester work on its own where it won't risk confusing
you. We assume REGTEST_SCRATCH is defined to point to this location.

2. Clone the required repos into the scratch area.

    ```
    git clone https://github.com/AMReX-Codes/amrex.git ${REGTEST_SCRATCH}/amrex
    git clone https://github.com/AMReX-Codes/regression_testing.git ${REGTEST_SCRATCH}/regression_testing
    git clone https://github.com/AMReX-Codes/incflo.git ${REGTEST_SCRATCH}/incflo
    ```

3.  Move to the location where the tests will be built/run, create
testing area expected by reg test scripts:

    ```
    cd ${REGTEST_SCRATCH}; mkdir test_data
    ```

4.  Edit the config file, ${REGTEST_SCRATCH}/incflo-tests.ini
to set the AMReX and incflo scratch clone locations and desired branch/SHA.
In that file, also set

     testTopDir =  ${REGTEST_SCRATCH}/test_data/incflo

Example config files are available in the regression_testing repo. Hopefully we will add an incflo one here soon. 

5. (optional) Create symbolic link
    ```
    ln -s ${REGTEST_SCRATCH}/regression_testing/regtest.py .
    ```

6.  Generate the initial benchmark solutions for all the tests listed
in the .ini configuration file.  Rerunning this at any time will
overwrite the previous versions of the benchmarks

    ```
    ./regtest.py --make_benchmarks "<a useful comment>" incflo-tests.ini
    ```

5. Upon some trigger event, re-run the tests and format the results in
html.  In this case, the results will appear as
test_data/incflo/web/index.html

    ```
    ./regtest.py incflo-tests.ini
    ```

NOTE: The regtest.py script takes a handy option "--no_update All",
which instructs the tester to work with the scratch repositories as
they currently exist.  Without this option specified, the branch of
each repository that is specified in the .ini config file is "git
pull"'d to obtain its most recent version; the original state of the
repositories are restored when the tests complete.  Using this
feature, a user can checkout any specific branch of any of the
repositories in the scratch area and run the complete set of tests.  A
user may wish to do this prior to issuing a "pull request", for
example.

More information on available options is given by
    ```
    ./regtest.py -h
    ```
which prints a verbose description of usage and setup. 
