target_sources(${amr_wind_lib_name} PRIVATE

  RefinementCriteria.cpp
  CartBoxRefinement.cpp
  FieldRefinement.cpp
  QCriterionRefinement.cpp
  VorticityRefinement.cpp
  OversetRefinement.cpp
  GeometryRefinement.cpp

  # Geometry based refinement options
  BoxRefiner.cpp
  CylinderRefiner.cpp
  )
