#ifndef _INCFLO_BCLIST_H_
#define _INCFLO_BCLIST_H_

// Boundary conditions list auxiliary structure
//
class BCList {
  public:
    BCList(const int undefined = 0,
           const int ig = 9,
           const int pinf = 10,
           const int pout = 11,
           const int minf = 20,
           const int nsw = 100,
           const int slip = 101,
           const int wall_model = 102)
      : und_(undefined)
      , ig_(ig)
      , pinf_(pinf)
      , pout_(pout)
      , minf_(minf)
      , nsw_(nsw)
      , slip_(slip)
      , wall_model_(wall_model)
    {}

    int get_undefined() const {return und_;};
    int get_ig() const {return ig_;};
    int get_pinf() const {return pinf_;};
    int get_pout() const {return pout_;};
    int get_minf() const {return minf_;};
    int get_nsw() const {return nsw_;};
    int get_slip() const {return slip_;};
    int get_wall_model() const {return wall_model_;};

  private:
    const int und_;
    const int ig_;
    const int pinf_;
    const int pout_;
    const int minf_;
    const int nsw_;
    const int slip_;
    const int wall_model_;
};

#endif
