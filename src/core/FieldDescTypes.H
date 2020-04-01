#ifndef FIELDDESCTYPES_H
#define FIELDDESCTYPES_H

namespace amr_wind {

enum class FieldState : int {
    NP1 = 0,
    N,
    NM1,
    NPH,
    NMH,
    New = NP1,
    Old = N,
};

enum class FieldLoc : int { CELL, NODE, XFACE, YFACE, ZFACE };

} // namespace amr_wind

#endif /* FIELDDESCTYPES_H */
