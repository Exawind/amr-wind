.. _compiling:

Compiling
==================

[Outdated documentation on compiling](https://exawind.github.io/amr-wind/user/build.html)

[Spack-Manager example for ExaWind Codes](https://sandialabs.github.io/spack-manager/user_profiles/developers/snapshot_workflow.html)

[Spack-Manager abbreviated example](https://sandialabs.github.io/spack-manager/user_profiles/developers/developer_workflow.html#quick-start)

Before you run AMR-Wind, you must first compile the LES solver. You may also need to compile additional codes in tandem (e.g., OpenFAST if you want to model turbines). There are two general approaches to compiling AMR-Wind: building from source or using [Spack-Manager](https://github.com/psakievich/spack-manager) ([docs](https://sandialabs.github.io/spack-manager/index.html)). I recommend Spack-Manager, because it can get very complicated to manually compile with all the interdependencies.

In my mind, [Spack](https://github.com/spack/spack) and Spack-Manager give us something like Python's `conda install <library>`, except for software packages. Spack is a hefty and powerful tool, and Spack-Manager was developed to help simplify it for ExaWind users.

If you're trying to quickly set up AMR-Wind, the [Snapshot Developer Workflow Example](https://sandialabs.github.io/spack-manager/user_profiles/developers/snapshot_workflow.html) is a quick way to cover the basics. I'll present a streamlined process here, but if things break, check out the official docs.

First, download Spack-Manager. I recommend setting it up somewhere where you have lots of storage (e.g., on Eagle, a `/projects/` directory instead of `/home/`), because it's easy to run out of space quickly. Spack-Manager also makes it much easier to install AMR-Wind on your own machine, if you are interested in running small cases locally. If you use a Mac, the compiler option "apple-clang" is recommended, instead of "gcc" as is used in this walkthrough.

```
export SPACKMANDIR=<fill in the directory>
cd ${SPACKMANDIR}
git clone --recursive https://github.com/sandialabs/spack-manager.git
```

Then activate it:
```
export SPACK_MANAGER=${SPACKMANDIR}/spack-manager && source ${SPACK_MANAGER}/start.sh && spack-start
```

Before we install anything, let's get some info about the codes that Spack can see
```
spack info amr-wind
```
and
```
spack info openfast
```
With these commands, you can see info like the safe versions as well as different flags that can get specified.

Also, before we compile, let's load the version of `gcc` that we will be using. As an example, let's suppose that `gcc/11.2.0` is the compiler version available on our machine via a module.
```
module load gcc/11.2.0
```

With this done, let's now create our environment. Environments, which are a grouping of specifications for a particular installation of the code, can be made to "belong" to a directory, using the `-d` option, or belong to a name, using the `-n` option. Environments created with a name will have their files copied into `$SPACK_MANAGER/environments`. In my experience, environments belonging to a directory are easier to modify and access.
```
mkdir ${SPACKMANDIR}/spack-july2023
cd ${SPACKMANDIR}/spack-july2023
quick-create-dev -d ${SPACKMANDIR}/spack-july2023 -s openfast@master%gcc@11.2.0 amr-wind@main+openfast+hdf5%gcc@11.2.0
```
This `quick-create-dev` command has flags selected so that that AMR-Wind will work with OpenFAST, and AMR-Wind also has the option to save out certain files using HDF5. If you forget the `+openfast` flag, your AMR-Wind simulations of turbines will crash and give you a confusing error message.

As part of the creation of the environment, Spack-Manager will automatically clone the repositories listed (amr-wind and openfast). If you already have a clone that you would like to use (or your own fork), you can instead create a soft/symbolic link within the directory prior to the environment creation step, and no new cloning will take place. Just ensure that the link has the exact name of the library as Spack-Manager understands it ("amr-wind" or "openfast").

The repositories cloned by Spack-Manager are shallow clones, and do not automatically have any commit history. If you would like to compile an older version of a code using a different commit, you can retrieve the commit history using the command `git fetch --unshallow` within the repository and then check out any past commit that you may need. After choosing a different commit, be sure to run `git submodule update` to modify the submodules to correspond to the chosen commit.

Note: the fact we specified "master" and "main" branches when we created the environment does not mean that the code in these repositories must be on the master and main branches, respectively. These references communicate to Spack-Manager a grouping of dependencies for each code. In almost all cases, using different commits or even your own fork for these codes does not change their dependencies, and so the specification of "master" and "main" is the correct spec for whatever version of the code you are using.

With the source code present and in the state you desire, let Spack compile everything by running
```
spack install
```
This process will take some time, and it will generate a large wall of text. Once that command is done running, we can also confirm the completion of the process by locating the amr-wind executable. Within the amr-wind directory, there will be a folder named `spack-build-<hash>`, where the <hash> is a signature of letters and numbers assigned to this repository and this environment. Within this build directory will be the executable `amr_wind`. 


### Additional details for turbine simulations
If you are simulating an OpenFAST turbine, you will probably need a controller for your turbine. I recommend [ROSCO](https://github.com/NREL/ROSCO). The full installation docs are [here](https://rosco.readthedocs.io/en/latest/source/install.html#full-rosco-installation), but we only need the controller from ROSCO, so the installation process is light. Below, I'll document what I do on Eagle.

0. `module load` the libraries that you will be loading when you run AMR-Wind, e.g. on Eagle `module load gcc/8.4.0 mpt cmake` (note: I run with `gcc/8.4.0` instead of `11.2.0` like I demoed for DelftBlue)
1. `cd` into a directory (ideally somewhere in `/projects/`), and clone ROSCO `git clone https://github.com/NREL/ROSCO.git`
2. Compile ROSCO

```
# Compile ROSCO
cd ROSCO/ROSCO
mkdir build
cd build
cmake ..                        # Mac/linux only
make install
```

This will generate a file called `libdiscon.so`, and we will later point to this file in OpenFAST configuration files.