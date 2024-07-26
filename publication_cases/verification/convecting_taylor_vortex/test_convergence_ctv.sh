$mpiexec_hybrid -np 1 ~/amr-wind/spack-build-vscvwnv/amr_wind ctv_godunov.inp amr.n_cell=   8   8 4 geometry.prob_hi = 2. 2. 0.5    geometry.prob_lo = 0. 0. -0.5    > output_8.txt
tail -1 ctv.log
mv ctv.log ctv_8.log
$mpiexec_hybrid -np 1 ~/amr-wind/spack-build-vscvwnv/amr_wind ctv_godunov.inp amr.n_cell=  16  16 4 geometry.prob_hi = 2. 2. 0.25   geometry.prob_lo = 0. 0. -0.25   > output_16.txt
tail -1 ctv.log
mv ctv.log ctv_16.log
$mpiexec_hybrid -np 4 ~/amr-wind/spack-build-vscvwnv/amr_wind ctv_godunov.inp amr.n_cell=  32  32 4 amr.blocking_factor_x = 8 amr.blocking_factor_y = 8 geometry.prob_hi = 2. 2. 0.125  geometry.prob_lo = 0. 0. -0.125  > output_32.txt
tail -1 ctv.log
mv ctv.log ctv_32.log
$mpiexec_hybrid -np 4 ~/amr-wind/spack-build-vscvwnv/amr_wind ctv_godunov.inp amr.n_cell=  64  64 4 amr.blocking_factor_x = 16 amr.blocking_factor_y = 16 geometry.prob_hi = 2. 2. 0.0625 geometry.prob_lo = 0. 0. -0.0625 > output_64.txt
tail -1 ctv.log
mv ctv.log ctv_64.log
$mpiexec_hybrid -np 8 ~/amr-wind/spack-build-vscvwnv/amr_wind ctv_godunov.inp amr.n_cell= 128 128 4 amr.blocking_factor_x = 16 amr.blocking_factor_y = 16 geometry.prob_hi = 2. 2. 0.03125 geometry.prob_lo = 0. 0. -0.03125 > output_128.txt
tail -1 ctv.log
mv ctv.log ctv_128.log
