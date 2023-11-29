/**This is the header file that should be used when compiling a standalone
 * shared library for an ABLCodedInlet.
 *
 * A minimal example (ConstantInlet.cpp):

#include "coded_inlet_fields.h"

void velocity(double time, double x, double y, double z, double (&ptvel)[3]) { ptvel[0]=10; ptvel[1]=0; ptvel[2]=0; }
void temperature(double time, double x, double y, double z, double &pttemp) { pttemp = 290; }

 * Compiled (on a mac) with:
 *
 *   clang++ -v ConstantInlet.cpp -o ConstantInlet.dylib -shared -fPIC
 *
 * To activate, set in the input file:
 *
 *   codedInlet.lib = "ConstantInlet.dylib"
 *
 */

extern "C" {

void velocity(double time, double x, double y, double z, double (&ptvel)[3]);
void temperature(double time, double x, double y, double z, double &pttemp);

}
