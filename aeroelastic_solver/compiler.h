#pragma once

#if defined(__clang__)

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
#define VECTORIZE_LOOP _Pragma("clang loop vectorize_width(2)")
#else
#define VECTORIZE_LOOP _Pragma("clang loop vectorize_width(scalable)")
#endif

#elif defined(__GNUC__)
// GCC does not support vectorize_width pragma, so use ivdep
#define VECTORIZE_LOOP _Pragma("GCC ivdep")

#else
#define VECTORIZE_LOOP

#endif
