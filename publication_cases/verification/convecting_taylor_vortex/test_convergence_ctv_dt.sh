$mpiexec_hybrid -np 4 ~/amr-wind/spack-build-vscvwnv/amr_wind ctv_godunov.inp amr.n_cell=  32  32 4 geometry.prob_hi = 2. 2. 0.125  geometry.prob_lo = 0. 0. -0.125  time.cfl = 0.05 time.stop_time = 3.0 time.max_step  = 100 > output_32_CFL005.txt
tail -1 ctv.log
mv ctv.log ctv_32_CFL005.log

$mpiexec_hybrid -np 4 ~/amr-wind/spack-build-vscvwnv/amr_wind ctv_godunov.inp amr.n_cell=  32  32 4 geometry.prob_hi = 2. 2. 0.125  geometry.prob_lo = 0. 0. -0.125  time.cfl = 0.1 time.stop_time = 3.0 time.max_step  = 100 > output_32_CFL01.txt
tail -1 ctv.log
mv ctv.log ctv_32_CFL01.log

$mpiexec_hybrid -np 4 ~/amr-wind/spack-build-vscvwnv/amr_wind ctv_godunov.inp amr.n_cell=  32  32 4 geometry.prob_hi = 2. 2. 0.125  geometry.prob_lo = 0. 0. -0.125  time.cfl = 0.25 time.stop_time = 3.0 time.max_step  = 100 > output_32_CFL025.txt
tail -1 ctv.log
mv ctv.log ctv_32_CFL025.log

$mpiexec_hybrid -np 4 ~/amr-wind/spack-build-vscvwnv/amr_wind ctv_godunov.inp amr.n_cell=  32  32 4 geometry.prob_hi = 2. 2. 0.125  geometry.prob_lo = 0. 0. -0.125  time.cfl = 0.5 time.stop_time = 3.0 time.max_step  = 100 > output_32_CFL05.txt
tail -1 ctv.log
mv ctv.log ctv_32_CFL05.log

$mpiexec_hybrid -np 4 ~/amr-wind/spack-build-vscvwnv/amr_wind ctv_godunov.inp amr.n_cell=  32  32 4 geometry.prob_hi = 2. 2. 0.125  geometry.prob_lo = 0. 0. -0.125  time.cfl = 0.75 time.stop_time = 3.0 time.max_step  = 100 > output_32_CFL075.txt
tail -1 ctv.log
mv ctv.log ctv_32_CFL075.log

$mpiexec_hybrid -np 4 ~/amr-wind/spack-build-vscvwnv/amr_wind ctv_godunov.inp amr.n_cell=  32  32 4 geometry.prob_hi = 2. 2. 0.125  geometry.prob_lo = 0. 0. -0.125  time.cfl = 0.95 time.stop_time = 3.0 time.max_step  = 100 > output_32_CFL095.txt
tail -1 ctv.log
mv ctv.log ctv_32_CFL095.log



