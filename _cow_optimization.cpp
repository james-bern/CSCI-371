struct HessianEntry {
    int row;
    int col;
    double val;
};

double *sparse2dense(int R, int C, StretchyBuffer<HessianEntry> *sparse) {
    double *dense = (double *) calloc(R * C, sizeof(double));
    #define RXC(M, row, col) ((M)[C * (row) + (col)])
    for_(k, sparse->length) { RXC(dense, sparse->data[k].row, sparse->data[k].col) += sparse->data[k].val; }
    #undef RXC
    return dense;
}
#ifdef EIGEN_LINEAR_SOLVER
struct EigenTriplet { int row, col; real val; };

void eigenSimplicialCholesky(int N, double *x, int A_length, EigenTriplet *A_data, double *b);
#endif
// x = inv(A) * b, where A is NxN
void solve_sparse_linear_system(int N, double *x, StretchyBuffer<HessianEntry> *_A, double *b) {
    { // checks
        ASSERT(x);
        ASSERT(N);
        ASSERT(_A);
        ASSERT(b);
    }

    if (_A->length == 0) {
        printf("A empty\n");
        memset(x, 0, N * sizeof(double));
        return;
    }

    #ifdef USE_EIGEN
    {
        eigenSimplicialCholesky(N, x, _A->length, (EigenTriplet *) _A->data, b);
    }
    #else
    {
        do_once { printf("[warn] USE_EIGEN not #define'd; falling back to dense gauss-jordan\n"); };

        // build the augmented matrix
        double *A = sparse2dense(N, N + 1, _A);
        #define NXNP1(M, row, col) ((M)[(N + 1) * (row) + (col)])
        for_(k, N) { NXNP1(A, k, N) = b[k]; }
        defer { free(A); };

        { // convert to triangular form (in place)
            // https://en.wikipedia.org/wiki/Gaussian_elimination
            int m = N;
            int n = N + 1;
            int h = 0;
            int k = 0;

            double *scratch = (double *) malloc(n * sizeof(double));
            defer { free(scratch); };

            while (h < m && k < n) {
                int max_i = -1;
                double max_abs = -INFINITY;
                {
                    for (int i = h; i < m; ++i) {
                        double tmp = ABS(NXNP1(A, i, k));
                        if (tmp > max_abs) {
                            max_abs = tmp;
                            max_i = i;
                        }
                    }
                }
                if (IS_EQUAL(0., NXNP1(A, max_i, k))) {
                    ++k;
                } else {
                    { // for_(c, n) { SWAP(NXNP1(A, h, c), NXNP1(A, max_i, c)); }
                        double *row_a = A + n * h;
                        double *row_b = A + n * max_i;
                        int size = n * sizeof(double);
                        memcpy(scratch, row_a, size);
                        memcpy(row_a, row_b, size);
                        memcpy(row_b, scratch, size);
                    }
                    for (int i = h + 1; i < m; ++i) {
                        double f = NXNP1(A, i, k) / NXNP1(A, h, k);
                        NXNP1(A, i, k) = 0;
                        for (int j = k + 1; j < n; ++j) {
                            NXNP1(A, i, j) = NXNP1(A, i, j) - NXNP1(A, h, j) * f;
                        }
                    }
                    ++h;
                    ++k;
                }
            }
        }

        // back substitue and store result in x
        {
            memset(x, 0, N * sizeof(double));
            for (int row = N - 1; row >= 0; --row) {
                for (int col = N - 1; col >= row; --col) {
                    x[row] += NXNP1(A, row, col) * NXNP1(A, col, N);
                }
            }
        }
        #undef NXNP1
    }
    #endif
}

double Vector_dot(int N, double *u, double *v) {
    double ret = 0;
    for_(i, N) {
        ret += u[i] * v[i];
    }
    return ret;
}

// convenience functions to add/extract segments and blocks
inline void add(double *a, int i, vec2 a_i) {
    a[2 * i + 0] += a_i.x;
    a[2 * i + 1] += a_i.y;
};
inline void add(StretchyBuffer<HessianEntry> *A, int i, int j, mat2 A_ij) {
    sbuff_push_back(A, { 2 * i + 0, 2 * j + 0, A_ij(0, 0) });
    sbuff_push_back(A, { 2 * i + 1, 2 * j + 0, A_ij(1, 0) });
    sbuff_push_back(A, { 2 * i + 1, 2 * j + 1, A_ij(1, 1) });
    sbuff_push_back(A, { 2 * i + 0, 2 * j + 1, A_ij(0, 1) });
};

