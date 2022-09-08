// msvc details switches /Bt /d2cgsummary 

#define real double
#define for_(i, N) for (int i = 0; i < N; ++i)

#ifndef NELEMS
#define NELEMS(fixed_size_array) int(sizeof(fixed_size_array) / sizeof((fixed_size_array)[0]))
#endif
#ifndef ASSERT
#define ASSERT(b) do { if (!(b)) { printf("ASSERT Line %d in %s\n", __LINE__, __FILE__); printf("press Enter to crash"); getchar(); *((volatile int *) 0) = 0; } } while (0)
#endif

#define NUM_DENm1(f, F) (double(f) / ((F) - 1))
#define NUM_DEN(f, F) (double(f) / (F))

#define SWAP(a, b) do {                   \
    ASSERT(sizeof(a) == sizeof(b));       \
    void *__SWAP_tmp = malloc(sizeof(a)); \
    memcpy(__SWAP_tmp, &a, sizeof(a));    \
    memcpy(&a, &b, sizeof(b));            \
    memcpy(&b, __SWAP_tmp, sizeof(b));    \
} while (0)

// // http://nothings.org/stb_ds/
// *((int*)(stb)-8)
#define STB_DS_IMPLEMENTATION
#include "ext/stb/stb_ds.h"
#undef arrput
#undef arrlen
#undef arrlenu
template <typename T> void arrput(T *&array, T element) { stbds_arrput(array, element); }
#define arrlen(array) int(stbds_arrlen(array))
template <typename T> T *jim_stb2raw(T *dyn) {
    int size = arrlen(dyn) * sizeof(T);
    T *raw = (T *) malloc(size);
    memcpy(raw, dyn, size);
    arrfree(dyn);
    return raw;
}
template <typename T> T *jim_raw2stb(int n, T *raw) {
    T *dyn = 0;
    for_(i, n) arrput(dyn, raw[i]);
    return dyn;
}

// https://handmade.network/forums/t/1273-post_your_c_c++_macro_tricks/3
#define __defer(line) defer_ ## line
#define _defer(line) __defer(line)
#define defer auto _defer(__LINE__) = defer_dummy() + [&]( )
template <typename F> struct Defer { Defer(F f) : f(f) {} ~Defer() { f(); } F f; }; template <typename F> Defer<F> makeDefer( F f ) { return Defer<F>( f ); }; struct defer_dummy {}; template<typename F> Defer<F> operator+( defer_dummy, F&& f ) { return makeDefer<F>( std::forward<F>(f) ); }


// https://en.cppreference.com/w/c/algorithm/qsort
void jim_sort_against(void *base, int nitems, int size, real *corresp_values_to_sort_against) {
    struct qsortHelperStruct {
        int index;
        real value;
    };
    qsortHelperStruct *helperArray = (qsortHelperStruct *) calloc(sizeof(qsortHelperStruct), nitems); {
        for_(i, nitems) helperArray[i] = { i, corresp_values_to_sort_against[i] };
    }

    int(* comp)(qsortHelperStruct *, qsortHelperStruct *) \
        = [](qsortHelperStruct *a, qsortHelperStruct *b) -> int { return (a->value < b->value) ? -1 : 1; }; 

    qsort(helperArray, nitems, sizeof(qsortHelperStruct), (int (*)(const void *, const void *))comp);

    void *tmp_buffer = malloc(nitems * size); { // fornow
        for_(i, nitems) memcpy(\
                ((char *) tmp_buffer) + (i * size), \
                ((char *) base) + (helperArray[i].index * size), \
                size);
        memcpy(base, tmp_buffer, nitems * size);
    } free(tmp_buffer);
}
