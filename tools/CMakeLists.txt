set ( FCOMPARE_PAR_EXE fcompare_par )
set ( FJOIN_PAR_EXE    fjoin_par )
set ( FEXTRACT_EXE     fextract )

add_executable( ${FCOMPARE_PAR_EXE} EXCLUDE_FROM_ALL F90/fcompare_par.f90 )

add_executable( ${FJOIN_PAR_EXE} EXCLUDE_FROM_ALL F90/fjoin_par.f90 )

add_executable( ${FEXTRACT_EXE} EXCLUDE_FROM_ALL C++/fextract.cpp )
target_link_libraries( ${FEXTRACT_EXE} AMReX::amrex ) 

# Communicate to main Cmakelists file
set ( INCFLO_F90_TOOLS ${FEXTRACT_EXE} ${FEXTRACT_EXE} ${FJOIN_PAR_EXE} PARENT_SCOPE  )
set ( FEXTRACT_EXE ${FEXTRACT_EXE} PARENT_SCOPE )
set ( FCOMPARE_PAR_EXE ${FCOMPARE_PAR_EXE} PARENT_SCOPE )
set ( FJOIN_PAR_EXE ${FJOIN_PAR_EXE} PARENT_SCOPE )
