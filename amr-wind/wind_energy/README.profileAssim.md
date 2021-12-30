# Profile Assimilation Development Notes
Eliot Quon, Shashank Yellapantula

TODOs:
- [ ] Implement general input for weighting matrix
- [ ] Generalize polynomial fit (i.e., implement `m_ind_polyOrder`)
- [ ] Generalize "Wrf" naming in code to "Meso"
- [ ] Deal with assumption that the ABL statistics class computes statistics at the cell-centeres
      only on level 0


## Differences from SOWFA

- Regression matrix (Z^T W Z) is formed with z/zmax in SOWFA
  Regression matrix (Z^T W Z) is formed with z*scaleFact in AMR-Wind, with hard-coded scaleFact=1e-3
- Blending/forcing transition in SOWFA is over:
    (assimMaxHeight, assimMaxHeight+blendThickness)
  Blending/forcing transition in AMR-Wind is over:
    (transition_height, transition_height+transition_thickness)

