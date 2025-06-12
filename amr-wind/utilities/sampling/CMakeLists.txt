target_sources(${amr_wind_lib_name}
  PRIVATE

  Sampling.cpp
  SamplingContainer.cpp
  SamplingUtils.cpp
  LineSampler.cpp
  LidarSampler.cpp
  RadarSampler.cpp
  DTUSpinnerSampler.cpp
  PlaneSampler.cpp
  ProbeSampler.cpp
  FreeSurfaceSampler.cpp
  VolumeSampler.cpp
  # The ones below should probably be moved
  FieldNorms.cpp
  KineticEnergy.cpp
  Enstrophy.cpp
  WaveEnergy.cpp
  )
