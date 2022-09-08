#define ASSERT(b) do { if (!(b)) {                        \
    printf("ASSERT Line %d in %s\n", __LINE__, __FILE__); \
    printf("press Enter to crash"); getchar();                \
    *((volatile int *) 0) = 0;                                \
} } while (0)
#define STATIC_ASSERT(cond) static_assert(cond, "STATIC_ASSERT");

#define PI 3.14159265358
#define RAD(deg) (PI / 180 * (deg))
#define DEG(rad) (180. / PI * (rad))

#define ABS(a) ((a) < 0 ? -(a) : (a))
#define SGN(a) ((a) < 0 ? -1 : 1)
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define AVG(a, b) (.5 * (a) + .5 * (b))
#define CLAMP(t, a, b) MIN(MAX(t, a), b)

#define LERP(t, a, b) ((1 - (t)) * (a) + (t) * (b))
#define INVERSE_LERP(p, a, b) (((p) - (a)) / double((b) - (a)))

#define TINY 1e-7
#define IS_POSITIVE(a) ((a) > TINY)
#define IS_NEGATIVE(a) ((a) < -TINY)
#define IS_ZERO(a) (ABS(a) < TINY)
#define IS_EQUAL(a, b) (ABS((a) - (b)) < TINY)
#define IN_RANGE(c, a, b) (((a) - TINY < (c)) && ((c) < (b) + TINY))

#define IS_EVEN(a) ((a) % 2 == 0)
#define IS_ODD(a) ((a) % 2 != 0)

#define MODULO(x, N) (((x) % (N) + (N)) % (N)) // works on negative numbers
#define NELEMS(fixed_size_array) int(sizeof(fixed_size_array) / sizeof((fixed_size_array)[0]))

#define STR(foo) #foo
#define XSTR(foo) STR(foo)
#define CONCAT_(a, b) a ## b
#define CONCAT(a, b) CONCAT_(a, b)

// do_once (this is an unorthodox macro, but we need it to work around a vs code bug)
#define do_once \
    static bool CONCAT(_do_once_, __LINE__) = false; \
    bool CONCAT(_prev_do_once_, __LINE__) = CONCAT(_do_once_, __LINE__); \
    CONCAT(_do_once_, __LINE__) = true; \
    if (!CONCAT(_prev_do_once_, __LINE__) && CONCAT(_do_once_, __LINE__))

#define INCHES(mm) ((mm) / 25.4)
#define MM(inches) ((inches) * 25.4)
