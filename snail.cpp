// { snail.h }                           .----.   @   @                         
//                                      / .-"-.`.  \v/                          
//                                      | | '\ \ \_/ )                          
//  a smol and slow linalg library    ,-\ `-.' /.'  /                           
//  ---------------------------------'---`----'----'            james-bern 2022 

// vec2, vec3, and vec4 are   2-,   3-, and   4-vectors  respectively           
// mat2, mat3, and mat4 are 2x2-, 3x3-, and 4x4-matrices respectively           
//                                                                              
// all arithmetic operators you might want are defined                          
// e.g. for vec3's s, v and mat3 M, the line s += -(M * v) / 2; is valid        
//                                                                              
// V2(x, y), V3(x, y, z), and V4(x, y, z, w) are "constructors"                 
// vec3 v = { 1, 2, 3 };     <- don't need V3 if just creating a vec3           
// vec3 v = 2 * V3(1, 2, 3); <-  _do_ need V3 if creating and using on same line
// when in doubt, just use V2, V3, V4                                           
//                                                                              
// the data of vector v can be accessed as                                      
// -  v[0], v[1], v[2], v[3]                                                    
// -   v.x,  v.y,  v.y,  v.w                                                    
// -   v.r,  v.g,  v.b,  v.a                                                    
// -  v.xy                                                                      
// - v.xyz                                                                      
// - or if you prefer, v.data[0], v.data[1], ...                                
//                                                                              
// the data of TxT matrix M can be accessed as M(r, c)                          
// - or if you prefer, M.data[T * r + c]                                        
//                                                                              
// snail_linalg.h has no explicit notion of column vectors vs. row vectors      
// -   M v is coded M * v                                                       
// - v^T M is coded v * M                                                       
//                                                                              
// the * operator is _not_ overloaded for two vectors                           
// - inner product a^T b is coded   dot(a, b)                                   
// - outer product a b^T is coded outer(a, b)                                   
//                                                                              
// inverse(M) returns the inverse of M, determinant(M) returns determinant...   

#ifndef SNAIL_WAS_INCLUDED
#define SNAIL_WAS_INCLUDED

// dependencies ////////////////////////////////////////////////////////////////

#define _CRT_SECURE_NO_WARNINGS
#include <cstdio>
#include <cmath>
#include <cstring>

// internal macros /////////////////////////////////////////////////////////////

#define SNAIL_ASSERT(b) do { if (!(b)) {                        \
    printf("ASSERT Line %d in %s\n", __LINE__, __FILE__); \
    printf("press Enter to crash"); getchar();                \
    *((volatile int *) 0) = 0;                                \
} } while (0)
#define SNAIL_FOR_(i, N) for (int i = 0; i < N; ++i)

// vectors and matrices ////////////////////////////////////////////////////////

template <int T> union SnailVec {
    double data[T];
    double &operator [](int index) { return data[index]; }
};

template <int T> union SnailMat {
    double data[T * T];
    double &operator ()(int r, int c) { return data[T * r + c]; }
};

// sugary accessors ////////////////////////////////////////////////////////////

template <> union SnailVec<2> {
    struct { double x, y; };
    double data[2];
    double &operator [](int index) { return data[index]; }
};
template <> union SnailVec<3> {
    struct { double x, y, z; };
    struct { double r, g, b; };
    struct { SnailVec<2> xy; double _; };
    double data[3];
    double &operator [](int index) { return data[index]; }
};
template <> union SnailVec<4> {
    struct { double x, y, z, w; };
    struct { double r, g, b, a; };
    struct { SnailVec<3> xyz; double _; };
    double data[4];
    double &operator [](int index) { return data[index]; }
};

// "constructors" //////////////////////////////////////////////////////////////

SnailVec<2> snail_V2(double x, double y) { SnailVec<2> ret = { x, y }; return ret; }
SnailVec<3> snail_V3(double x, double y, double z) { SnailVec<3> ret = { x, y, z }; return ret; }
SnailVec<4> snail_V4(double x, double y, double z, double w) { SnailVec<4> ret = { x, y, z, w }; return ret; }

// snort names /////////////////////////////////////////////////////////////////

#define V2 snail_V2
#define V3 snail_V3
#define V4 snail_V4

#define Vec SnailVec
#define Mat SnailMat

typedef SnailVec<2> vec2;
typedef SnailVec<3> vec3;
typedef SnailVec<4> vec4;
typedef SnailMat<2> mat2;
typedef SnailMat<3> mat3;
typedef SnailMat<4> mat4;

// arithmetic operators ////////////////////////////////////////////////////////

// vectors
template <int T> SnailVec<T>  operator +  (SnailVec<T> A, SnailVec<T> B) {
    SnailVec<T> result;
    SNAIL_FOR_(i, T) {
        result.data[i] = A.data[i] + B.data[i];
    }
    return result;
}
template <int T> SnailVec<T> &operator += (SnailVec<T> &A, SnailVec<T> B) {
    A = A + B;
    return A;
}

template <int T> SnailVec<T>  operator -  (SnailVec<T> A, SnailVec<T> B) {
    SnailVec<T> result;
    SNAIL_FOR_(i, T) {
        result.data[i] = A.data[i] - B.data[i];
    }
    return result;
}
template <int T> SnailVec<T> &operator -= (SnailVec<T> &A, SnailVec<T> B) {
    A = A - B;
    return A;
}

template <int T> SnailVec<T>  operator *  (double scalar, SnailVec<T> v) {
    SnailVec<T> result;
    SNAIL_FOR_(i, T) {
        result.data[i]  = scalar * v.data[i];
    }
    return result;
}
template <int T> SnailVec<T>  operator *  (SnailVec<T> v, double scalar) {
    SnailVec<T> result = scalar * v;
    return result;
}
template <int T> SnailVec<T> &operator *= (SnailVec<T> &v, double scalar) {
    v = scalar * v;
    return v;
}
template <int T> SnailVec<T>  operator -  (SnailVec<T> v) {
    return -1 * v;
}

template <int T> SnailVec<T>  operator /  (SnailVec<T> v, double scalar) {
    SnailVec<T> result;
    SNAIL_FOR_(i, T) {
        result.data[i]  = v.data[i] / scalar;
    }
    return result;
}
template <int T> SnailVec<T> &operator /= (SnailVec<T> &v, double scalar) {
    v = v / scalar;
    return v;
}

// matrices
template <int T> SnailMat<T>  operator +  (SnailMat<T> A, SnailMat<T> B) {
    SnailMat<T> ret = {};
    SNAIL_FOR_(k, T * T) {
        ret.data[k] = A.data[k] + B.data[k];
    }
    return ret;
}
template <int T> SnailMat<T> &operator += (SnailMat<T> &A, SnailMat<T> B) {
    A = A + B;
    return A;
}

template <int T> SnailMat<T>  operator -  (SnailMat<T> A, SnailMat<T> B) {
    SnailMat<T> ret = {};
    SNAIL_FOR_(i, T * T) {
        ret.data[i] = A.data[i] - B.data[i];
    }
    return ret;
}
template <int T> SnailMat<T> &operator -= (SnailMat<T> &A, SnailMat<T> B) {
    A = A + B;
    return A;
}

template <int T> SnailMat<T>  operator *  (SnailMat<T> A, SnailMat<T> B) {
    SnailMat<T> ret = {};
    SNAIL_FOR_(r, T) {
        SNAIL_FOR_(c, T) {
            SNAIL_FOR_(i, T) {
                ret(r, c) += A(r, i) * B(i, c);
            }
        }
    }
    return ret;
}
template <int T> SnailMat<T> &operator *= (SnailMat<T> &A, SnailMat<T> B) {
    A = A * B;
    return A;
}
template <int T> SnailVec<T>  operator *  (SnailMat<T> A, SnailVec<T> b) { // A b
    SnailVec<T> ret = {};
    SNAIL_FOR_(r, T) {
        SNAIL_FOR_(c, T) {
            ret[r] += A(r, c) * b[c];
        }
    }
    return ret;
}
template <int T> SnailVec<T>  operator *  (SnailVec<T> b, SnailMat<T> A) { // b^T A
    SnailVec<T> ret = {};
    SNAIL_FOR_(r, T) {
        SNAIL_FOR_(c, T) {
            ret[r] += A(c, r) * b[c];
        }
    }
    return ret;
}
template <int T> SnailMat<T>  operator *  (double scalar, SnailMat<T> M) {
    SnailMat<T> result = {};
    SNAIL_FOR_(k, T * T) {
        result.data[k] = scalar * M.data[k];
    }
    return result;
}
template <int T> SnailMat<T>  operator *  (SnailMat<T> M, double scalar) {
    return scalar * M;
}
template <int T> SnailMat<T> &operator *= (SnailMat<T> &M, double scalar) {
    M = scalar * M;
    return M;
}
template <int T> SnailMat<T>  operator -  (SnailMat<T> M) {
    return -1 * M;
}

template <int T> SnailMat<T>  operator /  (SnailMat<T> M, double scalar) {
    return (1 / scalar) * M;
}
template <int T> SnailMat<T> &operator /= (SnailMat<T> &M, double scalar) {
    M = M / scalar;
    return M;
}

// important vector functions //////////////////////////////////////////////////

template <int T> double dot(SnailVec<T> A, SnailVec<T> B) {
    double result = 0;
    for (int i = 0; i < T; ++i) {
        result += A.data[i] * B.data[i];
    }
    return result;
}
template <int T> SnailMat<T> outer(SnailVec<T> u, SnailVec<T> v) {
    SnailMat<T> ret = {};
    SNAIL_FOR_(r, T) {
        SNAIL_FOR_(c, T) {
            ret(r, c) = u[r] * v[c];
        }
    }
    return ret;
}

double cross(SnailVec<2> A, SnailVec<2> B) {
    return A.x * B.y - A.y * B.x;
}
SnailVec<3> cross(SnailVec<3> A, SnailVec<3> B) {
    return { A.y * B.z - A.z * B.y, A.z * B.x - A.x * B.z, A.x * B.y - A.y * B.x };
}

template <int T> double squaredNorm(SnailVec<T> v) {
    return dot(v, v);
}
template <int T> double norm(SnailVec<T> v) {
    return sqrt(squaredNorm(v));
}
template <int T> SnailVec<T> normalized(SnailVec<T> v) {
    double norm_v = norm(v);
    SNAIL_ASSERT(fabs(norm_v) > 1e-7);
    return (1 / norm_v) * v;
}

// important matrix functions //////////////////////////////////////////////////

template <int T> SnailMat<T> transpose(SnailMat<T> M) {
    SnailMat<T> ret = {};
    SNAIL_FOR_(r, T) {
        SNAIL_FOR_(c, T) {
            ret(r, c) = M(c, r);
        }
    }
    return ret;
}

double determinant(SnailMat<2> M) {
    return M(0, 0) * M(1, 1) - M(0, 1) * M(1, 0);
}
double determinant(SnailMat<3> M) {
    return M(0, 0) * (M(1, 1) * M(2, 2) - M(2, 1) * M(1, 2))
        - M(0, 1) * (M(1, 0) * M(2, 2) - M(1, 2) * M(2, 0))
        + M(0, 2) * (M(1, 0) * M(2, 1) - M(1, 1) * M(2, 0));
}
double determinant(SnailMat<4> M) {
    double A2323 = M(2, 2) * M(3, 3) - M(2, 3) * M(3, 2);
    double A1323 = M(2, 1) * M(3, 3) - M(2, 3) * M(3, 1);
    double A1223 = M(2, 1) * M(3, 2) - M(2, 2) * M(3, 1);
    double A0323 = M(2, 0) * M(3, 3) - M(2, 3) * M(3, 0);
    double A0223 = M(2, 0) * M(3, 2) - M(2, 2) * M(3, 0);
    double A0123 = M(2, 0) * M(3, 1) - M(2, 1) * M(3, 0);
    return M(0, 0) * ( M(1, 1) * A2323 - M(1, 2) * A1323 + M(1, 3) * A1223 ) 
        - M(0, 1) * ( M(1, 0) * A2323 - M(1, 2) * A0323 + M(1, 3) * A0223 ) 
        + M(0, 2) * ( M(1, 0) * A1323 - M(1, 1) * A0323 + M(1, 3) * A0123 ) 
        - M(0, 3) * ( M(1, 0) * A1223 - M(1, 1) * A0223 + M(1, 2) * A0123 ) ;
}

SnailMat<2> inverse(SnailMat<2> M) {
    double invdet = 1 / determinant(M);
    return { invdet * M(1, 1), 
        invdet * -M(0, 1), 
        invdet * -M(1, 0), 
        invdet * M(0, 0) };
}
SnailMat<3> inverse(SnailMat<3> M) {
    double invdet = 1 / determinant(M);
    return { invdet * (M(1, 1) * M(2, 2) - M(2, 1) * M(1, 2)),
        invdet * (M(0, 2) * M(2, 1) - M(0, 1) * M(2, 2)),
        invdet * (M(0, 1) * M(1, 2) - M(0, 2) * M(1, 1)),
        invdet * (M(1, 2) * M(2, 0) - M(1, 0) * M(2, 2)),
        invdet * (M(0, 0) * M(2, 2) - M(0, 2) * M(2, 0)),
        invdet * (M(1, 0) * M(0, 2) - M(0, 0) * M(1, 2)),
        invdet * (M(1, 0) * M(2, 1) - M(2, 0) * M(1, 1)),
        invdet * (M(2, 0) * M(0, 1) - M(0, 0) * M(2, 1)),
        invdet * (M(0, 0) * M(1, 1) - M(1, 0) * M(0, 1)) };
}
SnailMat<4> inverse(SnailMat<4> M) {
    double invdet = 1 / determinant(M);
    double A2323 = M(2, 2) * M(3, 3) - M(2, 3) * M(3, 2) ;
    double A1323 = M(2, 1) * M(3, 3) - M(2, 3) * M(3, 1) ;
    double A1223 = M(2, 1) * M(3, 2) - M(2, 2) * M(3, 1) ;
    double A0323 = M(2, 0) * M(3, 3) - M(2, 3) * M(3, 0) ;
    double A0223 = M(2, 0) * M(3, 2) - M(2, 2) * M(3, 0) ;
    double A0123 = M(2, 0) * M(3, 1) - M(2, 1) * M(3, 0) ;
    double A2313 = M(1, 2) * M(3, 3) - M(1, 3) * M(3, 2) ;
    double A1313 = M(1, 1) * M(3, 3) - M(1, 3) * M(3, 1) ;
    double A1213 = M(1, 1) * M(3, 2) - M(1, 2) * M(3, 1) ;
    double A2312 = M(1, 2) * M(2, 3) - M(1, 3) * M(2, 2) ;
    double A1312 = M(1, 1) * M(2, 3) - M(1, 3) * M(2, 1) ;
    double A1212 = M(1, 1) * M(2, 2) - M(1, 2) * M(2, 1) ;
    double A0313 = M(1, 0) * M(3, 3) - M(1, 3) * M(3, 0) ;
    double A0213 = M(1, 0) * M(3, 2) - M(1, 2) * M(3, 0) ;
    double A0312 = M(1, 0) * M(2, 3) - M(1, 3) * M(2, 0) ;
    double A0212 = M(1, 0) * M(2, 2) - M(1, 2) * M(2, 0) ;
    double A0113 = M(1, 0) * M(3, 1) - M(1, 1) * M(3, 0) ;
    double A0112 = M(1, 0) * M(2, 1) - M(1, 1) * M(2, 0) ;
    return { invdet * ( M(1, 1) * A2323 - M(1, 2) * A1323 + M(1, 3) * A1223 ),
        invdet * - ( M(0, 1) * A2323 - M(0, 2) * A1323 + M(0, 3) * A1223 ),
        invdet *   ( M(0, 1) * A2313 - M(0, 2) * A1313 + M(0, 3) * A1213 ),
        invdet * - ( M(0, 1) * A2312 - M(0, 2) * A1312 + M(0, 3) * A1212 ),
        invdet * - ( M(1, 0) * A2323 - M(1, 2) * A0323 + M(1, 3) * A0223 ),
        invdet *   ( M(0, 0) * A2323 - M(0, 2) * A0323 + M(0, 3) * A0223 ),
        invdet * - ( M(0, 0) * A2313 - M(0, 2) * A0313 + M(0, 3) * A0213 ),
        invdet *   ( M(0, 0) * A2312 - M(0, 2) * A0312 + M(0, 3) * A0212 ),
        invdet *   ( M(1, 0) * A1323 - M(1, 1) * A0323 + M(1, 3) * A0123 ),
        invdet * - ( M(0, 0) * A1323 - M(0, 1) * A0323 + M(0, 3) * A0123 ),
        invdet *   ( M(0, 0) * A1313 - M(0, 1) * A0313 + M(0, 3) * A0113 ),
        invdet * - ( M(0, 0) * A1312 - M(0, 1) * A0312 + M(0, 3) * A0112 ),
        invdet * - ( M(1, 0) * A1223 - M(1, 1) * A0223 + M(1, 2) * A0123 ),
        invdet *   ( M(0, 0) * A1223 - M(0, 1) * A0223 + M(0, 2) * A0123 ),
        invdet * - ( M(0, 0) * A1213 - M(0, 1) * A0213 + M(0, 2) * A0113 ),
        invdet *   ( M(0, 0) * A1212 - M(0, 1) * A0212 + M(0, 2) * A0112 ) };
}

// using 4x4 transforms ////////////////////////////////////////////////////////

template <int T> SnailVec<T> transform_point(const SnailMat<4> &M, SnailVec<T> p) {
    SnailVec<4> p_hom = {};
    memcpy(p_hom.data, p.data, T * sizeof(double));
    p_hom.w = 1;
    SnailVec<4> ret_hom = M * p_hom;
    ret_hom /= ret_hom.w;
    SnailVec<T> ret = {};
    memcpy(ret.data, ret_hom.data, T * sizeof(double));
    return ret;
}
template <int T> SnailVec<T> transform_direction(const SnailMat<4> &M, SnailVec<T> p) {
    SnailVec<4> p_hom = {};
    memcpy(p_hom.data, p.data, T * sizeof(double));
    SnailVec<4> ret_hom = M * p_hom;
    SnailVec<T> ret = {};
    memcpy(ret.data, ret_hom.data, T * sizeof(double));
    return ret;
}
template <int T> SnailVec<T> transform_normal(const SnailMat<4> &M, SnailVec<T> p) {
    SnailVec<4> p_hom = {};
    memcpy(p_hom.data, p.data, T * sizeof(double));
    SnailVec<4> ret_hom = inverse(transpose(M)) * p_hom;
    SnailVec<T> ret = {};
    memcpy(ret.data, ret_hom.data, T * sizeof(double));
    return ret;
}

// 4x4 transform cookbook //////////////////////////////////////////////////////

template <int T> SnailMat<T> IdentityMatrix() {
    SnailMat<T> ret = {};
    for (int i = 0; i < T; ++i) {
        ret(i, i) = 1;
    }
    return ret;
}
const SnailMat<4> Identity = IdentityMatrix<4>();
SnailMat<4> Translation(double x, double y, double z = 0) {
    SnailMat<4> ret = Identity;
    ret(0, 3) = x;
    ret(1, 3) = y;
    ret(2, 3) = z;
    return ret;
}
SnailMat<4> Translation(SnailVec<2> xy) {
    return Translation(xy.x, xy.y);
}
SnailMat<4> Translation(SnailVec<3> xyz) {
    return Translation(xyz.x, xyz.y, xyz.z);
}
SnailMat<4> Scaling(double s) {
    SnailMat<4> ret = {};
    ret(0, 0) = s;
    ret(1, 1) = s;
    ret(2, 2) = s;
    ret(3, 3) = 1;
    return ret;
}
SnailMat<4> Scaling(double x, double y, double z = 1) {
    SnailMat<4> ret = {};
    ret(0, 0) = x;
    ret(1, 1) = y;
    ret(2, 2) = z;
    ret(3, 3) = 1;
    return ret;
}
SnailMat<4> Scaling(SnailVec<2> xy) {
    return Scaling(xy.x, xy.y);
}
SnailMat<4> Scaling(SnailVec<3> xyz) {
    return Scaling(xyz.x, xyz.y, xyz.z);
}
SnailMat<4> RotationX(double t) {
    SnailMat<4> ret = Identity;
    ret(1, 1) = cos(t); ret(1, 2) = -sin(t);
    ret(2, 1) = sin(t); ret(2, 2) =  cos(t);
    return ret;
}
SnailMat<4> RotationY(double t) {
    SnailMat<4> ret = Identity;
    ret(0, 0) =  cos(t); ret(0, 2) = sin(t);
    ret(2, 0) = -sin(t); ret(2, 2) = cos(t);
    return ret;
}
SnailMat<4> RotationZ(double t) {
    SnailMat<4> ret = Identity;
    ret(0, 0) = cos(t); ret(0, 1) = -sin(t);
    ret(1, 0) = sin(t); ret(1, 1) =  cos(t);
    return ret;
}
SnailMat<4> Rotation(SnailVec<3> axis, double angle) {
    double x = axis.x;
    double y = axis.y;
    double z = axis.z;
    double x2 = x * x;
    double y2 = y * y;
    double z2 = z * z;
    double xy = x * y;
    double xz = x * z;
    double yz = y * z;
    double c = cos(angle);
    double s = sin(angle);
    double d = 1-c;
    return { c+x2*d, xy*d-z*s, xz*d+y*s, 0,
        xy*d+z*s, c+y2*d, yz*d-x*s, 0,
        xz*d-y*s, yz*d+x*s, c+z2*d, 0,
        0, 0, 0, 1 };
}

// optimization stuff //////////////////////////////////////////////////////////

template <int T> SnailMat<T> firstDerivativeofUnitVector(SnailVec<T> v) {
    SnailVec<T> tmp = normalized(v);
    return (1 / norm(v)) * (IdentityMatrix<T>() - outer(tmp, tmp));
}
#define firstDerivativeOfNorm normalized
#define secondDerivativeOfNorm firstDerivativeofUnitVector

template <int T> double squaredNorm(SnailMat<T> M) {
    double ret = 0;
    for (int i = 0; i < T * T; ++i) {
        ret += M.data[i] * M.data[i];
    }
    return ret;
}

// misc functions //////////////////////////////////////////////////////////////

template <int T> SnailVec<T> cwiseAbs(SnailVec<T> A) {
    for (int i = 0; i < T; ++i) A[i] = abs(A[i]);
    return A;
}
SnailVec<2> e_theta(double theta) {
    return { cos(theta), sin(theta) };
}
SnailVec<2> perpendicularTo(SnailVec<2> v) {
    return { v.y, -v.x };
}

// utility /////////////////////////////////////////////////////////////////////

template <int T> void pprint(SnailVec<T> v) {
    printf("[ ");
    SNAIL_FOR_(i, T) {
        printf("%lf", v[i]);
        if (i != T - 1) printf(", ");
    }
    printf(" ]\n");
}
template <int T> void pprint(SnailMat<T> M) {
    SNAIL_FOR_(r, T) {
        printf("| ");
        SNAIL_FOR_(c, T) {
            printf("%lf", M(r, c));
            if (c != T - 1) printf(", ");
        }
        printf(" |\n");
    }
}

#undef SNAIL_FOR_
#undef SNAIL_ASSERT
#endif
