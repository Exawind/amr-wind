.. _getting-started:

Getting started with AMR-Wind
================================

Discretization and solution approach
------------------------------------

The numerical methdology used to solve the partial differential equations (PDEs)
within AMR-Wind is documented in `Almgren et al. (JCP 1998)
<https://ccse.lbl.gov/Publications/almgren/abchw.pdf>`_.

Navigating source code
------------------------

``AMR-Wind`` is built on top of `AMReX library
<https://amrex-codes.github.io/amrex/docs_html/>`_. Users are strongly
recommended to read through the AMReX documentation and understand the basic
AMReX concepts before jumping into the AMR-Wind source code.

The `Basics section
<https://amrex-codes.github.io/amrex/docs_html/Basics_Chapter.html>`_ provides a
thorough overview of the basic datastructures and ways to interact with these
structures. The `GPU section
<https://amrex-codes.github.io/amrex/docs_html/GPU_Chapter.html>`_ provides an
overview of the AMReX GPU strategy and the higher-level functions (e.g.,
``parallel-for`` abstractions) available to write GPU-ready code within
AMR-Wind. The `Linear Solvers section
<https://amrex-codes.github.io/amrex/docs_html/LinearSolvers_Chapter.html>`_
provides an overview of the multi-level multigrid (MLMG) solvers used to solve
the various linear systems within AMR-Wind.
