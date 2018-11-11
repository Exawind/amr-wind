#ifndef ALGOIM_REAL_HPP
#define ALGOIM_REAL_HPP

// Header file for defining Algoim::Real. In most use cases, Algoim::Real is a typedef for double. When
// high precisiona arithmetic is used via the QD library, Algoim::Real may instead be a double-double or
// quadruple-double number -- in such cases, a (global) compiler preprocessor definition for ALGOIM_DDREAL
// or ALGOIM_QDREAL must be defined to enable QD interoperability.

#if defined(ALGOIM_QDREAL)
#include <qd/qd_real.h>
#include <qd/fpu.h>
#define ALGOIM_HPREAL
#elif defined(ALGOIM_DDREAL)
#include <qd/dd_real.h>
#include <qd/fpu.h>
#define ALGOIM_HPREAL
#else
#include <cstdio>
#endif

// If using QD, inject some of its functions into the std namespace
#ifdef ALGOIM_HPREAL
namespace std
{
    inline Real abs(const Real& x) { return ::abs(x); }
    inline Real sqrt(const Real& x) { return ::sqrt(x); }
}
#endif

namespace Algoim
{
    // Typedefs for Algoim::Real
    #if defined(ALGOIM_QDREAL)
    typedef qd_real Real;
    #elif defined(ALGOIM_DDREAL)
    typedef dd_real Real;
    #else
    typedef double Real;
    #endif

    // Definition for strToReal, which takes a C-style char array and returns a Real
    #ifdef ALGOIM_HPREAL
    inline Real strToReal(const char* s)
    {
        return Real(s);
    }
    #else
    inline double strToReal(const char* s)
    {
        double x;
        std::sscanf(s, "%lf", &x);
        return x;
    }
    #endif
}

#endif
