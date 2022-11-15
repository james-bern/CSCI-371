#define _CRT_SECURE_NO_WARNINGS
#include "snail.cpp"
#include "cow.cpp"
#include "_cow_supplement.cpp"
#include "jim.cpp"
// #define EIGEN_LINEAR_SOLVER // will discuss Thursday
#include "_cow_optimization.cpp"

struct {
    bool playing = false;
    bool check_derivatives = false;
    int num_rows = 2;
    bool noisy = false;
    double theta_gravity = -RAD(90);
} tweaks = {};

struct Simulation {
    // num_*
    int num_nodes;
    int num_springs;
    int num_pins;

    // data
    double *X; // node_rest_positions
    int2 *springs;
    int *pins;

    // convenience
    int num_bytes_x; // sim.X, state.x, state.x_prev, state.x_prev_prev
    int N; // 2 * num_nodes
};

struct State {
    double *x;           // node_current_positions           (x_k  )
    double *x_prev;      // node_previous_positions          (x_km1)
    double *x_prev_prev; // node_previous_previous_positions (x_km2)
};

struct Parameters {
    double totalMass = 1.;
    double gravitationalConstant_NOTE_positive = 9.81;
    double springSpringConstant = 1e2;
    double pinSpringConstant = 1e3;
    double timestep = 1. / 60;
};

#define ENABLED_NUM_FLAGS 4
struct Enabled {
    // NOTE in practice, i would probably use a bit field here #237x371CrossoverEpisode
    union {
        struct {
            bool gravity;
            bool springs;
            bool pins;
            bool dynamics;
        };
        bool flags[ENABLED_NUM_FLAGS];
    };
};
char *enabled_flag_padded_strings[ENABLED_NUM_FLAGS] = {
    " gravity",
    " springs",
    "    pins",
    "dynamics"
};


void simulation_draw(mat4 PV, Simulation *sim, State *state, Parameters *, Enabled *enabled) {
    gl_PV(PV);

    int num_nodes = sim->num_nodes;
    int num_springs = sim->num_springs;
    int num_pins = sim->num_pins;
    vec2 *X = (vec2 *) sim->X;
    vec2 *x = (vec2 *) state->x;
    int2 *springs = sim->springs;
    int *pins = sim->pins;

    { // springs
        // TODO springs as LINES                             
        //      dark blue if enabled->springs; otherwise gray


    }
    { // nodes
        // TODO nodes as blue POINTS


    }
    { // pins
        // TODO pins as POINTS (at x[pin]) and LINES (from x[pin] to X[pin])
        //      yellow if enabled->pins; otherwise gray                     


    }

    FORNOW_UNUSED(enabled);
    FORNOW_UNUSED(num_nodes);
    FORNOW_UNUSED(num_springs);
    FORNOW_UNUSED(num_pins);
    FORNOW_UNUSED(X);
    FORNOW_UNUSED(x);
    FORNOW_UNUSED(springs);
    FORNOW_UNUSED(pins);
}


// compute_and_add_to the energy U, gradient (dUdx) U_x, and Hessian (d2Udx2) U_xx
// note that U is more of an "energy" because of how we handle dynamics
void compute_and_add_to(Simulation *sim, State *state, Parameters *params, Enabled *enabled, double *U, double *U_x, StretchyBuffer<HessianEntry> *U_xx) {
    // // convenience
    // vec2 pointers to data
    vec2 *X           = (vec2 *) sim->X;
    vec2 *x           = (vec2 *) state->x;
    vec2 *x_prev      = (vec2 *) state->x_prev;
    vec2 *x_prev_prev = (vec2 *) state->x_prev_prev;
    vec2 *a = (vec2 *) malloc(sim->num_bytes_x); {
        for_(i, sim->num_nodes) { a[i] = (x[i] - 2 * x_prev[i] + x_prev_prev[i]) / pow(params->timestep, 2); }
    }
    defer { free(a); };

    if (enabled->gravity) {
        int num_nodes = sim->num_nodes;
        // vec2 *x;
        double m = params->totalMass / sim->num_nodes;
        double g = params->gravitationalConstant_NOTE_positive;

        // TODO implement gravity


        FORNOW_UNUSED(num_nodes);
        FORNOW_UNUSED(m);
        FORNOW_UNUSED(g);
    }

    if (enabled->springs) {
        // FORNOW HACK make diagonal springs less stiff
        double R_HACK = norm(X[0] - X[1]);

        double k_spring = params->springSpringConstant;
        for_(spring_i, sim->num_springs) {
            int i = sim->springs[spring_i].i;
            int j = sim->springs[spring_i].j;
            vec2 v = x[i] - x[j];
            double r = norm(v);
            vec2 r_v = firstDerivativeOfNorm(v);
            mat2 r_vv = secondDerivativeOfNorm(v);
            double R = norm(X[i] - X[j]);
            double Delta = r - R;

            // FORNOW HACK make diagonal springs less stiff
            if (R > TINY + R_HACK) { k_spring /= 16; }

            if (U) {
                (*U) += k_spring * pow(Delta, 2) / 2;
            }
            if (U_x) {
                vec2 seg = k_spring * Delta * r_v;
                add(U_x, i,  seg);
                add(U_x, j, -seg);
            }
            if (U_xx) {
                mat2 blk = k_spring * outer(r_v, r_v) + k_spring * Delta * r_vv;
                add(U_xx, i, i,  blk);
                add(U_xx, i, j, -blk);
                add(U_xx, j, j,  blk);
                add(U_xx, j, i, -blk);
            }
        }
    }

    if (enabled->pins) {
        int num_pins = sim->num_pins;
        // vec2 *x;
        // vec2 *X;
        int *pins = sim->pins;
        double k_pin = params->pinSpringConstant;

        // TODO implement pins


        FORNOW_UNUSED(num_pins);
        FORNOW_UNUSED(pins);
        FORNOW_UNUSED(k_pin);
    }

    if (enabled->dynamics) {
        double h = params->timestep;
        double m = params->totalMass / sim->num_nodes;
        mat2 I = M2(1, 0, 0, 1);
        for_(i, sim->num_nodes) {
            if (U) { (*U) += pow(h, 2) / 2 * m * squaredNorm(a[i]); }
            if (U_x) { add(U_x, i,  m * a[i]); }
            if (U_xx) { add(U_xx, i, i, m / pow(h, 2) * I); }
        }
    } else {
        double b = .0001;
        mat2 I = M2(1, 0, 0, 1);
        for_(i, sim->num_nodes) {
            if (U_x) { add(U_x, i,  b * x[i]); }
            if (U_xx) { add(U_xx, i, i, b * I); }
        }
    }
}

void finite_difference_and_add_to(Simulation *sim, State *state, Parameters *params, Enabled *enabled, double *U_x, StretchyBuffer<HessianEntry> *U_xx, double fd_stepsize = 0) {
    if (IS_EQUAL(0., fd_stepsize)) {
        fd_stepsize = 1e-5;
    }

    if (U_x) {
        double left;
        double right;
        for_(k, sim->N) {
            double x_k_o = state->x[k]; {
                state->x[k] = x_k_o - fd_stepsize;
                left = 0;
                compute_and_add_to(sim, state, params, enabled, &left, NULL, NULL);

                state->x[k] = x_k_o + fd_stepsize;
                right = 0;
                compute_and_add_to(sim, state, params, enabled, &right, NULL, NULL);
            } state->x[k] = x_k_o;
            U_x[k] += (right - left) / (2 * fd_stepsize);
        }
    }

    if (U_xx) {
        double *left = (double *) malloc(sim->N * sizeof(double));
        double *right = (double *) malloc(sim->N * sizeof(double));
        defer {
            free(left);
            free(right);
        };
        for_(col, sim->N) {
            double x_c_0 = state->x[col]; {
                state->x[col] = x_c_0 - fd_stepsize;
                memset(left, 0, sim->N * sizeof(double));
                compute_and_add_to(sim, state, params, enabled, NULL, left, NULL);

                state->x[col] = x_c_0 + fd_stepsize;
                memset(right, 0, sim->N * sizeof(double));
                compute_and_add_to(sim, state, params, enabled, NULL, right, NULL);
            } state->x[col] = x_c_0;
            for_(row, sim->N) {
                double val = (right[row] - left[row]) / (2 * fd_stepsize);
                if (!IS_EQUAL(0., val)) sbuff_push_back(U_xx, { row, col, val } );
            }
        }
    }
}

void check_derivatives(Simulation *sim, State *state, Parameters *params, Enabled *enabled, double fd_stepsize = 0) {
    #define ABSOLUTE_ERROR_THRESHOLD .0001
    #define RELATIVE_ERROR_THRESHOLD .001
    #define TRIGGER_INVALID(a) (isnan(a) || isinf(a))
    #define TRIGGER_ERROR(error) (error > ABSOLUTE_ERROR_THRESHOLD && (2. * error / (ABS(a) + ABS(b))) > RELATIVE_ERROR_THRESHOLD)

    double *U_x    = (double *) calloc(sim->N, sizeof(double));
    double *U_x_fd = (double *) calloc(sim->N, sizeof(double));
    double *U_xx;
    double *U_xx_fd;
    {
        StretchyBuffer<HessianEntry> _U_xx = {};
        StretchyBuffer<HessianEntry> _U_xx_fd = {};
        {
            compute_and_add_to(sim, state, params, enabled, NULL, U_x, &_U_xx);
            finite_difference_and_add_to(sim, state, params, enabled, U_x_fd, &_U_xx_fd, fd_stepsize);
            U_xx = sparse2dense(sim->N, sim->N, &_U_xx);
            U_xx_fd = sparse2dense(sim->N, sim->N, &_U_xx_fd);
        }
        sbuff_free(&_U_xx);
        sbuff_free(&_U_xx_fd);
    }

    { // check
        { // U_x
            bool passes = true;
            for_(k, sim->N) {
                double a = U_x[k];
                double b = U_x_fd[k];
                double error = ABS(a - b);
                bool invalid = TRIGGER_INVALID(a) || TRIGGER_INVALID(b);
                bool wrong = TRIGGER_ERROR(error);
                if (invalid || wrong) {
                    if (passes) { printf(" -- U_x FAIL"); }
                    passes = false;
                    if (invalid) { printf("%3d: nan or inf\n", k); }
                    else { printf("%2d: | (%lf) - (%lf) | = %lf\n", k, a, b, error); }
                }
            }
            if (passes) { printf(" -- U_x PASS"); }
        }
        { // U_xx
            bool passes = true;
            for_(r, sim->N) for_(c, sim->N) {
                #define NXN(M, row, col) ((M)[(sim->N) * (row) + (col)])
                double a = NXN(U_xx, r, c);
                double b = NXN(U_xx_fd, r, c);
                #undef NXN
                double error = ABS(a - b);
                bool invalid = TRIGGER_INVALID(a) || TRIGGER_INVALID(b);
                bool wrong = TRIGGER_ERROR(error);
                if (invalid || wrong) {
                    if (passes) { printf(" -- U_xx FAIL"); }
                    passes = false;
                    if (invalid) { printf("%3d, %3d: nan or inf\n", r, c); }
                    else { printf("%2d, %2d: | (%lf) - (%lf) | = %lf\n", r, c, a, b, error); }
                }
            }
            if (passes) { printf(" -- U_xx PASS"); }
        }
        printf("\n");
    }

    free(U_x);
    free(U_x_fd);
    free(U_xx);
    free(U_xx_fd);
}

// integrate forward one timestep using Newton's method for minimization with line search
// U_xx searchDir = -U_x
void step(Simulation *sim, State *state, Parameters *params, Enabled *enabled) {
    // FORNOW scratch
    StretchyBuffer<HessianEntry> U_xx = {};
    double *minus_U_x         = (double *) malloc(sim->num_bytes_x);
    double *searchDirection   = (double *) malloc(sim->num_bytes_x);
    double *_next_x           = (double *) malloc(sim->num_bytes_x);
    double *_next_x_0_line_search = (double *) malloc(sim->num_bytes_x);
    defer {
        sbuff_free(&U_xx);
        free(minus_U_x);
        free(searchDirection);
        free(_next_x);
        free(_next_x_0_line_search);
    };

    // warm start optimization at current position
    memcpy(_next_x, state->x, sim->num_bytes_x);
    State next  = {}; {
        next.x           = _next_x;
        next.x_prev      = state->x;
        next.x_prev_prev = state->x_prev;
    }

    int iterationOfNewtonWithLineSearch = 0;
    while (true) {
        { // compute_and_add_to -U_x, U_xx
            U_xx.length = 0;
            memset(minus_U_x, 0, sim->num_bytes_x);
            compute_and_add_to(sim, &next, params, enabled, NULL, minus_U_x, &U_xx);
            for_(k, sim->N) { minus_U_x[k] *= -1; }
        }

        if (Vector_dot(sim->N, minus_U_x, minus_U_x) < .001) { // convergence check
            break;
        }

        if (iterationOfNewtonWithLineSearch++ > 50) {
            if (tweaks.noisy) { printf("phyiscs solve failed\n"); }
            break;
        }

        { // get searchDirection (and regularize as needed)
            // U_x(x + searchDirection) ~ U_x + U_xx searchDirection := 0
            // => searchDirection = solve { U_xx searchDirection = -U_x }
            int iterationofDynamicRegularization = 0;
            do {
                if (iterationofDynamicRegularization == 1) { if (tweaks.noisy) { printf("not a descent direction\n"); } }
                solve_sparse_linear_system(sim->N, (double *) searchDirection, &U_xx, minus_U_x);
                // memcpy(searchDirection, minus_U_x, sim->num_bytes_x);

                { // regularize Hessian (sloppily)
                    for_(k, sim->N) {
                        sbuff_push_back(&U_xx, { k, k, pow(10, -4 + int(iterationofDynamicRegularization)) });
                    }
                    ++iterationofDynamicRegularization;
                    if (iterationofDynamicRegularization == 20) {
                        if (tweaks.noisy) { printf("dynamic regularization failed\n"); }
                        break;
                    }
                }
            } while (Vector_dot(sim->N, searchDirection, minus_U_x) < 0);
        }

        { // line search
            double O_0 = 0;
            compute_and_add_to(sim, &next, params, enabled, &O_0, NULL, NULL);

            memcpy(_next_x_0_line_search, next.x, sim->num_bytes_x);

            int iterationOfLineSearch = 0;
            double stepSize = 1;
            while (1) {
                // x_next = x_next_0 + stepSize * searchDirection
                for_(k, sim->N) { next.x[k] = _next_x_0_line_search[k] + stepSize * searchDirection[k]; }

                double O_curr = 0;
                compute_and_add_to(sim, &next, params, enabled, &O_curr, NULL, NULL);

                if (O_curr < O_0 + TINY) { // line search succeeded
                    break;
                }

                if (++iterationOfLineSearch > 30) {
                    if (tweaks.noisy) { printf("line search failed.\n"); }
                    break;
                }

                stepSize /= 2;
            }
        }
    }

    // state <- next
    memcpy(state->x_prev_prev, state->x_prev, sim->num_bytes_x);
    memcpy(state->x_prev,      state->x,      sim->num_bytes_x);
    memcpy(state->x,           next.x,        sim->num_bytes_x);
}

Simulation build_beam(int num_rows, int num_cols) {
    Simulation sim = {};
    #define INDEX(row, col) ((row) * num_cols + (col))
    { // nodes (README)
        sim.num_nodes = num_rows * num_cols;
        sim.X = (double *) calloc(sim.num_nodes, sizeof(vec2));
        sim.num_bytes_x = sim.num_nodes * sizeof(vec2);
        sim.N = 2 * sim.num_nodes;
        {
            double S = 1. / (num_cols);
            for_(row, num_rows) {
                double y = row * S;
                for_(col, num_cols) {
                    double x = col * S;
                    ((vec2 *) sim.X)[INDEX(row, col)] = { x, y };
                }
            }
        }
    }
    { // springs
        StretchyBuffer<int2> _springs = {};
        {
            for_(row, num_rows) for_(col, num_cols - 1) sbuff_push_back(&_springs, { INDEX(row, col), INDEX(row, col + 1) });
            for_(row, num_rows - 1) for_(col, num_cols) sbuff_push_back(&_springs, { INDEX(row, col), INDEX(row + 1, col) });
            for_(row, num_rows - 1) for_(col, num_cols - 1) {
                sbuff_push_back(&_springs, { INDEX(row, col), INDEX(row + 1, col + 1) });
                sbuff_push_back(&_springs, { INDEX(row + 1, col), INDEX(row, col + 1) });
            }
        }
        sim.springs = _springs.data;
        sim.num_springs = _springs.length;
    }
    { // pins
        StretchyBuffer<int> _pins = {};
        {
            for_(j, num_rows) {
                sbuff_push_back(&_pins, INDEX(j, 0));
                sbuff_push_back(&_pins, INDEX(j, num_cols - 1));
            }
        }
        sim.pins = _pins.data;
        sim.num_pins = _pins.length;
    }
    #undef INDEX
    return sim;
}

void hw10() {
    init();

    Camera2D camera = { 2, 0, -.5 };

    Simulation sim = build_beam(tweaks.num_rows, 3 * tweaks.num_rows);
    State state = {}; {
        state.x           = (double *) malloc(sim.num_bytes_x);
        state.x_prev      = (double *) malloc(sim.num_bytes_x);
        state.x_prev_prev = (double *) malloc(sim.num_bytes_x);
        memcpy(state.x,           sim.X, sim.num_bytes_x);
        memcpy(state.x_prev,      sim.X, sim.num_bytes_x);
        memcpy(state.x_prev_prev, sim.X, sim.num_bytes_x);
    }
    Parameters params = {};
    Enabled enabled = {}; {
        for_(k, ENABLED_NUM_FLAGS) {
            enabled.flags[k] = true;
        }
    }

    while (begin_frame()) {
        camera_move(&camera);
        mat4 PV = camera_get_PV(&camera);
        gl_PV(PV);

        { // gui
            imgui_readout("num_nodes", &sim.num_nodes);
            imgui_readout("num_springs", &sim.num_springs);
            imgui_readout("num_pins", &sim.num_pins);
            imgui_checkbox("noisy", &tweaks.noisy);
            imgui_checkbox("check_derivatives", &tweaks.check_derivatives, 'a');
            imgui_checkbox("gravity", &enabled.gravity, 'g');
            imgui_checkbox("pins", &enabled.pins, 'f');
            imgui_checkbox("dynamics", &enabled.dynamics, 'd');
            imgui_checkbox("springs", &enabled.springs, 's');
            imgui_checkbox("playing", &tweaks.playing, 'p');
            imgui_slider("timestep", &params.timestep, .001, .1);
            imgui_slider("springSpringConstant", &params.springSpringConstant, 1e0, 1e3);
            imgui_slider("pinSpringConstant", &params.pinSpringConstant, 1e0, 1e3);
            imgui_slider("theta_gravity", &tweaks.theta_gravity, -RAD(90) -PI, -RAD(90) + PI, true);
            { // rebuild
                int _num_rows = tweaks.num_rows;
                imgui_slider("num_rows", &tweaks.num_rows, 1, 20, 'j', 'k');
                if (_num_rows != tweaks.num_rows) {
                    { // fornow
                        free(sim.X);
                        free(sim.springs);
                        free(sim.pins);
                    }
                    sim = build_beam(tweaks.num_rows, 3 * tweaks.num_rows);
                    { // fornow
                        state.x           = (double *) realloc(state.x,           sim.num_bytes_x);
                        state.x_prev      = (double *) realloc(state.x_prev,      sim.num_bytes_x);
                        state.x_prev_prev = (double *) realloc(state.x_prev_prev, sim.num_bytes_x);
                        memcpy(state.x,           sim.X, sim.num_bytes_x);
                        memcpy(state.x_prev,      sim.X, sim.num_bytes_x);
                        memcpy(state.x_prev_prev, sim.X, sim.num_bytes_x);
                    }
                }
            }
            { // reset
                if (imgui_button("reset", 'r')) {
                    memcpy(state.x,           sim.X, sim.num_bytes_x);
                    memcpy(state.x_prev,      sim.X, sim.num_bytes_x);
                    memcpy(state.x_prev_prev, sim.X, sim.num_bytes_x);
                }
            }
        }

        if (tweaks.playing) {
            if (tweaks.check_derivatives) {
                for_(k, ENABLED_NUM_FLAGS) {
                    Enabled tmp = {};
                    tmp.flags[k] = true;
                    printf("%s", enabled_flag_padded_strings[k]);
                    check_derivatives(&sim, &state, &params, &tmp);
                }
            }
            step(&sim, &state, &params, &enabled);
        }
        widget_drag(PV, sim.num_nodes, (vec2 *) state.x, 0, monokai.white);

        { // draw
            simulation_draw(PV, &sim, &state, &params, &enabled);

            { // gravity vector
                gl_color(monokai.yellow);
                gl_begin(LINES);
                vec2 s = { 1.3, 0 };
                vec2 t = s + .3 * e_theta(tweaks.theta_gravity);
                gl_vertex(s);
                gl_vertex(t);
                double eps = .05;
                gl_vertex(t);
                gl_vertex(t + eps * e_theta(tweaks.theta_gravity - RAD(145)));
                gl_vertex(t);
                gl_vertex(t + eps * e_theta(tweaks.theta_gravity - RAD(215)));
                gl_end();
            }
        }
    }
}


int main() {
    hw10();
    return 0;
}

