typedef double real;
#define GL_REAL GL_DOUBLE
#define _CRT_SECURE_NO_WARNINGS

#ifndef JIM_NO_SNAIL
#ifdef JIM_IS_JIM
#define SNAIL_I_SOLEMNLY_SWEAR_I_AM_UP_TO_NO_GOOD
#endif
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

/*
// (very) optional style guide; feel free to delete entirely
// however, if you're newer to C/C++ then this guide may help reduce bugs
#if !defined(JIM_IS_JIM) && !defined(COW_NO_STYLE_GUIDE)
#define static     error__cow_style__prefer__setup_loop_structure__over__static
#define new        error__cow_style__prefer__malloc_free__over__new_delete
#define delete     error__cow_style__prefer__malloc_free__over__new_delete
// #define malloc     error__cow_style__prefer_calloc_over_malloc
#define float      error__cow_style__prefer_real_over_float
#define double     error__cow_style__prefer_real_over_double
#define unsigned   error__cow_style__prefer__int__over__unsigned_int__prefer_u8_over_unsigned_char
#define unique_ptr error__cow_style__prefer_raw_pointer
#define weak_ptr   error__cow_style__prefer_raw_pointer
#define shared_ptr error__cow_style__prefer_raw_pointer
#endif
*/
