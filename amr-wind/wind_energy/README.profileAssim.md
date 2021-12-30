# Profile Assimilation Development Notes
Eliot Quon, Shashank Yellapantula

TODOs:
- [ ] Implement general input for weighting matrix
- [ ] Generalize polynomial fit (i.e., implement `m_ind_polyOrder`)


## Differences from SOWFA

- Regression matrix (Z^T W Z) is formed with z/zmax in SOWFA
  Regression matrix (Z^T W Z) is formed with z*scaleFact in AMR-Wind, with hard-coded scaleFact=1e-3

