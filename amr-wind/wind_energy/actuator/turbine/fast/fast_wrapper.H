#ifndef FAST_WRAPPER_H
#define FAST_WRAPPER_H

namespace exw_fast {
#ifdef AMR_WIND_USE_OPENFAST
extern "C" {
#include "FAST_Library.h"
}

extern void exw_fast_output_redirect(char*);

inline constexpr int fast_strlen() { return INTERFACE_STRING_LENGTH; }

#else

#define ErrID_None 0
#define ErrID_Info 1
#define ErrID_Warn 2
#define ErrID_Severe 3
#define ErrID_Fatal 4

struct OpFM_InputType
{};
struct OpFM_OutputType
{};

struct SC_DX_InputType
{};
struct SC_DX_OutputType
{};

inline constexpr int fast_strlen() { return 1025; }

inline void exw_fast_output_redirect(char* /*unused*/) {}

inline void
FAST_AllocateTurbines(int* /*unused*/, int* /*unused*/, char* /*unused*/)
{}
inline void FAST_DeallocateTurbines(int* /*unused*/, char* /*unused*/) {}
inline void
FAST_OpFM_Solution0(int* /*unused*/, int* /*unused*/, char* /*unused*/)
{}
inline void FAST_OpFM_Step(int* /*unused*/, int* /*unused*/, char* /*unused*/)
{}
inline void FAST_CreateCheckpoint(
    int* /*unused*/, char* /*unused*/, int* /*unused*/, char* /*unused*/)
{}

inline void FAST_HubPosition(
    int* /*unused*/,
    float* /*unused*/,
    float* /*unused*/,
    double* /*unused*/,
    int* /*unused*/,
    char* /*unused*/)
{}

inline void FAST_OpFM_Init(
    int* /*unused*/,
    double* /*unused*/,
    const char* /*unused*/,
    int* /*unused*/,
    int* /*unused*/,
    int* /*unused*/,
    int* /*unused*/,
    float* /*unused*/,
    float* /*unused*/,
    int* /*unused*/,
    int* /*unused*/,
    float* /*unused*/,
    int* /*unused*/,
    double* /*unused*/,
    int* /*unused*/,
    int* /*unused*/,
    int* /*unused*/,
    OpFM_InputType* /*unused*/,
    OpFM_OutputType* /*unused*/,
    SC_DX_InputType* /*unused*/,
    SC_DX_OutputType* /*unused*/,
    int* /*unused*/,
    char* /*unused*/)
{}

inline void FAST_OpFM_Restart(
    int* /*unused*/,
    char* /*unused*/,
    int* /*unused*/,
    double* /*unused*/,
    int* /*unused*/,
    int* /*unused*/,
    int* /*unused*/,
    OpFM_InputType* /*unused*/,
    OpFM_OutputType* /*unused*/,
    SC_DX_InputType* /*unused*/,
    SC_DX_OutputType* /*unused*/,
    int* /*unused*/,
    char* /*unused*/)
{}
#endif
} // namespace exw_fast

#endif /* FAST_WRAPPER_H */
