# Profile Assimilation Development Notes
Eliot Quon, Shashank Yellapantula

TODOs:
- [x] Implement general input for weighting matrix
- [x] Update z normaliation in regression matrix
- [ ] Generalize polynomial fit (i.e., implement `m_ind_polyOrder`)
- [ ] Generalize "Wrf" naming in code to "Meso"
- [ ] Deal with assumption that the ABL statistics class computes statistics at the cell-centers
      only on level 0
- [ ] Separate ABLWrfForcing class into a separate source file
- [ ] Add r-test


## Differences from SOWFA

- Regression matrix (Z^T W Z) is formed with z/zmax in SOWFA
  Regression matrix (Z^T W Z) is formed with z*scaleFact in AMR-Wind, which was previously
    hard-coded with scaleFact=1e-3
- Blending/forcing transition in SOWFA is over:
    (assimMaxHeight, assimMaxHeight+blendThickness)
  Blending/forcing transition in AMR-Wind is over:
    (transition_height, transition_height+transition_thickness)

