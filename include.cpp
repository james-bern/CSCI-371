#define _CRT_SECURE_NO_WARNINGS
typedef double real;
#define GL_REAL GL_DOUBLE
#ifndef JIM_NO_SNAIL
#include "codebase/snail.cpp"
#endif
#include "codebase/cow.cpp"
#ifdef JIM_IS_JIM
#include "codebase/jim.cpp"
#endif

#include <iostream>
#include <fenv.h>
#include <stdlib.h>
#include <stdarg.h>
#include <utility>
#include <cstdio>
#include <cstring>
#include <cmath>

#if 0
// // optional style guide (consider disabling after hw0 if you already enjoy modern C++)
#define new        error__cow_style__prefer__malloc_free__over__new_delete
#define delete     error__cow_style__prefer__malloc_free__over__new_delete
#define malloc     error__cow_style__prefer_calloc_over_malloc
#define float      error__cow_style__prefer_real_over_float
#define double     error__cow_style__prefer_real_over_double
#define unsigned   error__cow_style__prefer__int__over__unsigned_int__prefer_u8_over_unsigned_char
#define unique_ptr error__cow_style__prefer_raw_pointer
#define weak_ptr   error__cow_style__prefer_raw_pointer
#define shared_ptr error__cow_style__prefer_raw_pointer
#endif
