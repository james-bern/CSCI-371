// #define COW_PATCH_FRAMERATE
// #define COW_PATCH_FRAMERATE_SLEEP
#include "include.cpp"

////////////////////////////////////////////////////////////////////////////////
// pendulum simulation /////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#define PENDULUM_SIMULATION_MODE_EXPLICIT 0
#define PENDULUM_SIMULATION_MODE_SEMI_IMPLICIT 1
#define PENDULUM_SIMULATION_MODE_IMPLICIT 2
#define _PENDULUM_SIMULATION_MODE_COUNT 3
char *_pendulum_simulation_mode_names[] = {
    "explicit euler",
    "semi-implicit euler",
    "implicit_euler",
};


// // parameters
// h is the simulation timestep
// g is the (signed) gravitational acceleration
// L is the length

// // state is { theta, omega }
// theta is the angle of the pendulum, measured counter-clockwise from the negative y-axix
// omega is the angular velocity of the pendulum

// NOTE: alpha is the angular acceleration of the pendulum


// simulate the pendulum forward one timestep and return the result (the next state)
// curr state is \xi_k = { \theta_k, \omega_k }
// next state is \xi_{ k + 1 } = { \theta_{k + 1}, \omega_{k + 1} }
vec2 step_forward_in_time(vec2 curr, int mode, real h, real g, real L) {
    // alpha is the angular acceleration of the pendulum
    // this is a lambda (local) function you can call just like a regular function
    auto get_alpha = [&](real theta) {
        return g / L * sin(theta);
    };

    vec2 next = curr;
    if (mode == PENDULUM_SIMULATION_MODE_EXPLICIT) {
        // TODO: explicit Euler

    } else if (mode == PENDULUM_SIMULATION_MODE_SEMI_IMPLICIT) {
        // TODO: semi-implicit Euler

    } else if (mode == PENDULUM_SIMULATION_MODE_IMPLICIT) {
        // TODO: implicit Euler

    }
    return next;
}


void hw10a() {
    vec2 state = { RAD(90.0), 0.0 };
    int mode = 0;

    real h = 0.0167; 
    real g = -9.81;  
    real L = 1.0;    

    bool paused = false;
    Camera2D camera = { 5.35, 0, -.75 };
    while (cow_begin_frame()) {
        camera_move(&camera);
        mat4 PV = camera_get_PV(&camera);
        bool reset = false;

        if (!paused) {
            state = step_forward_in_time(state, mode, h, g, L);
        }

        { // gui
            gui_slider("h", &h, .001, .1);
            gui_slider("L", &L, .1, 5);
            gui_slider("g", &g, -24, 24);
            gui_printf("");
                gui_slider("mode", &mode, 0, _PENDULUM_SIMULATION_MODE_COUNT - 1, 'j', 'k', true);
                gui_printf("mode: %s", _pendulum_simulation_mode_names[mode]);
            gui_printf("");
            gui_checkbox("paused", &paused, 'p');
            gui_printf("");
            { // reset
                static vec2 _state = state;
                if (gui_button("reset", 'r')) {
                    reset = true;
                    state = _state;
                    // sbuff_free(&trace_positions);
                }
            }
            gui_printf("");
            gui_readout("theta", &state[0]);
            gui_readout("omega", &state[1]);
        }

        { // draw / drag
            vec2 p = L * V2(sin(state[0]), -cos(state[0]));

            { // trace
                #define TRACE_LENGTH 128
                static vec3 trace_colors[TRACE_LENGTH];
                static StretchyBuffer<vec2> trace_positions = {};
                if (reset) { sbuff_free(&trace_positions); }

                static bool initialized;
                if (!initialized) {
                    initialized = true;

                    for (int i = 0; i < TRACE_LENGTH; ++i) {
                        trace_colors[i] = color_plasma(real(i) / real(TRACE_LENGTH - 1));
                    }
                }

                if (trace_positions.length < TRACE_LENGTH) {
                    sbuff_push_back(&trace_positions, p);
                } else {
                    memmove(trace_positions.data, trace_positions.data + 1, (TRACE_LENGTH - 1) * sizeof(vec2));
                    trace_positions.data[TRACE_LENGTH - 1] = p;
                }
                soup_draw(PV, SOUP_LINE_STRIP, trace_positions.length, trace_positions.data, trace_colors);
            }

            { // pendulum
                eso_color(monokai.white);
                eso_begin(PV, SOUP_LINES);
                eso_vertex(V2(0, 0));
                eso_vertex(p);
                eso_end();

                eso_color(monokai.yellow);
                eso_begin(PV, SOUP_POINTS, 16.0);
                eso_vertex(p);
                eso_end();

                if (widget_drag(PV, 1, &p)) {
                    p = L * normalized(p);
                    state[1] = 0;
                    state[0] = RAD(90) + atan2(p.y, p.x);
                }
            }

            { // energy plot
                static Plot plot = { 128 };
                if (plot.num_traces == 0) {
                    plot_add_trace(&plot, 0.0, 32.0, monokai.blue);
                    plot_add_trace(&plot, 0.0, 32.0, monokai.red);
                    plot_add_trace(&plot, 0.0, 32.0, monokai.purple);
                }
                if (reset) {
                    plot_clear(&plot);
                }

                mat4 PV_plot = M4_Translation(.3, -.3) * M4_Scaling(.5, .5 * _window_get_aspect());
                if (!paused) {
                    real m = 1.0;
                    real PE = m * g * (cos(state[0]) / L - L);
                    real KE = m * pow(L * state[1], 2) / 2.0;
                    plot_data_point(&plot, 0, PE);
                    plot_data_point(&plot, 1, KE);
                    plot_data_point(&plot, 2, PE + KE);
                }
                plot_draw(&plot, PV_plot);
                text_draw(PV_plot, "time", { 1.05, 0.0 });
                text_draw(PV_plot, "energy", { 0.05, 1.0 });
            }
        }
    }
}



////////////////////////////////////////////////////////////////////////////////
// skeletal animation //////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#define SKELETON_ANIMATION_MODE_ONLY_ZEROTH 0
#define SKELETON_ANIMATION_MODE_ALL_EQUAL 1
#define SKELETON_ANIMATION_MODE_RIGID 2
#define SKELETAL_ANIMATION_MODE_SMOOTH 3
#define _SKELETAL_ANIMATION_MODE_COUNT 4
char *skeletal_animation_mode_names[] = {
    "only-zeroth",
    "all-equal",
    "rigid",
    "smooth",
};

#define NUM_BONES 4
#define NUM_NODES (32 * 2)

// // parameters
// weights[i][j] is the skinning weight of bone j on node i
// L[j] is the length of bone j

// // bind pose of the system (aka rest pose)
// // this is the state of the system we use to calculate the weights
// b_bind[j] is the bind position of the near end of bone j in the skeleton (in world coordinates)
// NOTE: b_bind[0] is at the world origin
// NOTE: b_bind[NUM_BONES] is at the far end of the final bone
// NOTE: we assume the bind pose of the skeleton is lying along the x-axis
//       so we do NOT define theta_bind (theta_bind is understood to be 0.0 for all bones)
// s_bind[i] is the bind position of node i in the skin (in world coordinates)

// b_bind[0]      b_bind[1]   b_bind[2]        b_bind[3]     b_bind[4]          
// |              |           |                |             |                  
// v              v           v                v             v                  
// o--------------o-----------o----------------o-------------o   -> world x-axis
// ^                                                                            
// |                                                                            
// world origin                                                                 

// // current state of the system
// b[j] is the current position of the near end of bone j in the skeleton (in world coordinates)
// NOTE: b[0] should stay at the world origin forever
// NOTE: b[NUM_BONES] is the position of the far end of the final bone
// theta[j] is the current angle of bone j in the skeleton (in world coordinates)
// s[i] is the current position of node i in the skin (in world coordinates)


void update_skeleton(vec2 *b, real *L, real *theta) {
    // TODO: calculate b
    // (after you get this right, you should see the skeleton)

}


void update_skin(vec2 *s, real **weights, vec2 *b, real *theta, vec2 *b_bind, vec2 *s_bind) {
    // TODO: calculate s
    // (after you get this right, you should see the skin)

}


// helper funcion for initialize_weights
// returns the minimum distance from point p to line segment [a, b]
real point_segment_distance(vec2 p, vec2 a, vec2 b) {
    real L2 = squaredNorm(a - b);
    if (IS_ZERO(L2)) { return norm(p - a); }
    return norm(p - LERP(MAX(0, MIN(1, dot(p - a, b - a) / L2)), a, b));
}

void initialize_weights(real **weights, int mode, vec2 *s_bind, vec2 *b_bind) {
    for (int i = 0; i < NUM_NODES; ++i) {
        real w_i[NUM_BONES] = {}; {
            // TODO: calculate the i-th node's skinning weights
            // NOTE: they start out cleared to zero
            // NOTE: don't worry about ensuring sum(w_i) is 1 (i do it for you below)

            if (mode == SKELETON_ANIMATION_MODE_ONLY_ZEROTH) {
                // bind the node completely to the 0-th bone
                w_i[0] = 1.0;
            }
            else if (mode == SKELETON_ANIMATION_MODE_ALL_EQUAL) {
                // all bones equally-weighted
                for (int j = 0; j < NUM_BONES; ++j) {
                    w_i[j] = 1.0;
                }
            }
            else if (mode == SKELETON_ANIMATION_MODE_RIGID) {
                // TODO: bind the node completely to its closest bone
                // (after you get this right, the skin should look "like a robot" / "rigid")

            }
            else if (mode == SKELETAL_ANIMATION_MODE_SMOOTH) {
                // TODO (tricky): smoothly blend the weights between nearby bones
                // (after you get this right, the skin should look "squishy" / "smooth")

            }
        }

        { // normalize weights (ensure they sum(w_i) is 1)
            real sum = 0;
            for (int j = 0; j < NUM_BONES; ++j) {
                sum += w_i[j];
            }

            for (int j = 0; j < NUM_BONES; ++j) {
                w_i[j] /= sum;
            }
        }

        { // write w_i -> weights[i]
            for (int j = 0; j < NUM_BONES; ++j) {
                weights[i][j] = w_i[j];
            }
        }
    }
}


void hw10b() {
    real **weights[_SKELETAL_ANIMATION_MODE_COUNT] = {}; {
        for (int mode = 0; mode < _SKELETAL_ANIMATION_MODE_COUNT; ++mode) {
            weights[mode] = (real **) calloc(NUM_NODES, sizeof(real *));
            for (int i = 0; i < NUM_NODES; ++i) {
                weights[mode][i] = (real *) calloc(NUM_BONES, sizeof(real));
            }
        }
    }
    real L[NUM_BONES] = {};
    real _phi[NUM_BONES] = {};
    real theta[NUM_BONES] = {};
    vec2 b[NUM_BONES + 1] = {};
    vec2 b_bind[NUM_BONES + 1] = {};
    vec2 s[NUM_NODES] = {};
    vec2 s_bind[NUM_NODES] = {};

    { // initialize
        // bones
        for (int j = 0; j < NUM_BONES; ++j) {
            L[j] = LERP(real(j) / MAX(1, NUM_BONES - 1), 2.5, 1);
            b_bind[j + 1] = b_bind[j] + V2(L[j], 0);
        }

        // skin
        {
            // bind pose as a loop
            // e.g., for 10 nodes:
            // 9 8 7 6 5          
            // 0 1 2 3 4          
            int k = 0;
            for (int sign = -1; sign <= 1; sign += 2) {
                for (int i = 0; i < NUM_NODES / 2; ++i) {
                    real f = real(i) / real(NUM_NODES / 2 - 1);
                    if (sign > 0) { f = 1 - f; }
                    real skeleton_total_length = b_bind[NUM_BONES].x;
                    s_bind[k++] = V2(skeleton_total_length * f, .3 * sign);
                }
            }
        }

        for (int mode = 0; mode < _SKELETAL_ANIMATION_MODE_COUNT; ++mode) {
            initialize_weights(weights[mode], mode, s_bind, b_bind);
        }
    }

    // tweaks
    bool draw_character = false;
    bool draw_bind_pose = false;
    bool hide_skeleton = false;
    bool hide_nodes = false;
    bool hide_plots = false;

    // debug play
    bool debug_play = false;
    real debug_time = 0;

    // keyframing system
    #define MAX_KEYFRAMES 1024
    real keyframes[MAX_KEYFRAMES][NUM_BONES] = {};
    int num_keyframes = 0;
    bool tween_keyframes = false;
    const int FRAMES_PER_KEYFRAME = 32;
    int frame = 0;

    Camera2D camera = { 10, 0, 2 };
    int mode = 0;
    while (cow_begin_frame()) {
        camera_move(&camera);
        mat4 PV = camera_get_PV(&camera);

        { // compute absolute angles
            for (int j = 0; j < NUM_BONES; ++j) {
                theta[j] = ((j == 0) ? 0.0 : theta[j - 1]) + _phi[j];
            }
        }

        update_skeleton(b, L, theta);
        update_skin(s, weights[mode], b, theta, b_bind, s_bind);

        { // keyframing
            // animate the skeleton (i.e., set _phi using whatever method you like)
            // METHOD 1: press 'p' to play a sinusoidal trajectory                                 
            // METHOD 2: press 's' to save keyframes; press 't' to play them back with lerp        
            if (tween_keyframes && num_keyframes > 0) {
                int frame_A = (frame / FRAMES_PER_KEYFRAME) % num_keyframes;
                int frame_B = (frame_A + 1) % num_keyframes;
                real t = real(frame % FRAMES_PER_KEYFRAME) / FRAMES_PER_KEYFRAME;
                for (int j = 0; j < NUM_BONES; ++j) {
                    _phi[j] = LERP(
                            t,
                            keyframes[frame_A][j],
                            keyframes[frame_B][j]);
                }
                ++frame;
            } else if (debug_play) {
                for (int j = 0; j < NUM_BONES; ++j) {
                    _phi[j] += .02 * cos((j + 1) * debug_time / 2);
                }
                debug_time += .0167;
            }
        }

        { // gui
            for (int j = 0; j < NUM_BONES; ++j) {
                char buffer[24];
                sprintf(buffer, "joint angle %d", j);
                gui_slider(buffer, _phi + j, -PI, PI, true);
            }
            gui_printf("");
            gui_slider("mode", &mode, 0, _SKELETAL_ANIMATION_MODE_COUNT - 1, 'j', 'k', true);
            gui_printf("mode: %s", skeletal_animation_mode_names[mode]);
            if (!hide_plots) {
                mat4 PV_plot = M4_Translation(0, .25) * M4_Scaling(.78, .5 * .5 * _window_get_aspect());
                static vec2 trace[NUM_NODES / 2];
                {
                    vec2 axes[] = { { 1, 0 }, { 0, 0 }, { 0, 1 } };
                    soup_draw(PV_plot, SOUP_LINE_STRIP, 3, axes, NULL, monokai.gray);
                    text_draw(PV_plot, "x_rest", axes[0]);
                    text_draw(PV_plot, "weight", axes[2], V3(1, 1, 1), 0, { 0, -24 });
                }
                for (int j = 0; j < NUM_BONES; ++j) {
                    for (int i = 0; i < NUM_NODES / 2; ++i) {
                        trace[i] = { real(i) / real(NUM_NODES / 2 - 1), weights[mode][i][j] };
                    }
                    soup_draw(PV_plot, SOUP_LINE_STRIP, NUM_NODES / 2, trace, NULL, color_kelly(j));
                }
            }
            gui_printf("");
            gui_checkbox("debug_play", &debug_play, 'p');
            gui_printf("");
            if (gui_button("save keyframe", 's') && num_keyframes < MAX_KEYFRAMES) {
                for (int j = 0; j < NUM_BONES; ++j) {
                    keyframes[num_keyframes][j] = _phi[j];
                }
                ++num_keyframes;
            }
            gui_checkbox("tween keyframes", &tween_keyframes, 't');
            gui_readout("num_keyframes", &num_keyframes);
            gui_printf("");
            { // reset
                static real _bone_relative_angles_0[NUM_BONES];
                do_once { memcpy(_bone_relative_angles_0, _phi, sizeof(_phi)); };
                if (gui_button("reset", 'r')) {
                    memcpy(_phi, _bone_relative_angles_0, sizeof(_phi));

                    debug_time = 0;

                    memset(keyframes, 0, sizeof(keyframes));
                    num_keyframes = 0;
                    tween_keyframes = false;
                    frame = 0;
                }
            }
            gui_printf("");
            gui_checkbox("draw_character", &draw_character, 'z');
            gui_checkbox("draw_bind_pose", &draw_bind_pose, 'x');
            gui_checkbox("hide_skeleton", &hide_skeleton, 'c');
            gui_checkbox("hide_nodes", &hide_nodes, 'v');
            gui_checkbox("hide_plots", &hide_plots, 'b');
        }

        { // draw
            if (!hide_skeleton) {
                eso_begin(PV, GL_LINES); {
                    for (int j = 0; j < NUM_BONES; ++j) {
                        eso_color(color_kelly(j));
                        eso_vertex(b[j]);
                        eso_vertex(b[j + 1]);
                    }
                } eso_end();
                soup_draw(PV, SOUP_POINTS, NUM_BONES + 1, b, NULL, monokai.white);
            }

            if (draw_bind_pose) {
                soup_draw(PV, SOUP_POINTS, NUM_BONES + 1, b_bind, NULL, monokai.blue);
                soup_draw(PV, SOUP_LINE_STRIP, NUM_BONES + 1, b_bind, NULL, monokai.blue);
                soup_draw(PV, SOUP_LINE_LOOP, NUM_NODES, s_bind, NULL, .5 * monokai.blue);
                soup_draw(PV, SOUP_POINTS, NUM_NODES, s_bind, NULL, monokai.blue, 8.0);
            }

            if (!hide_nodes) {
                static vec3 *node_colors = (vec3 *) malloc(NUM_NODES * sizeof(vec3));
                for (int i = 0; i < NUM_NODES; ++i) {
                    node_colors[i] = {};
                    for (int j = 0; j < NUM_BONES; ++j) {
                        node_colors[i] += weights[mode][i][j] * color_kelly(j);
                    }
                }
                soup_draw(PV, SOUP_LINE_LOOP, NUM_NODES, s, node_colors, {}, 4.0);
                soup_draw(PV, SOUP_POINTS, NUM_NODES, s, node_colors, {}, 9.0);
            }

            if (draw_character) {
                const int NUM_QUADS = (NUM_NODES / 2 - 1);
                static vec2 vertex_positions[4 * NUM_QUADS];
                static vec3 vertex_colors[4 * NUM_QUADS];
                for (int i = 0; i < NUM_QUADS; ++i) {
                    vertex_positions[4 * i + 0] = s[i];
                    vertex_positions[4 * i + 1] = s[i + 1];
                    vertex_positions[4 * i + 2] = s[NUM_NODES - 1 - (i + 1)];
                    vertex_positions[4 * i + 3] = s[NUM_NODES - 1 - (i)];
                    real f = real(i) / (NUM_QUADS - 1);
                    vec3 color = color_rainbow_swirl(f);
                    vertex_colors[4 * i + 0] = color;
                    vertex_colors[4 * i + 1] = color;
                    vertex_colors[4 * i + 2] = color;
                    vertex_colors[4 * i + 3] = color;
                }
                soup_draw(PV, SOUP_QUADS, 4 * NUM_QUADS, vertex_positions, vertex_colors);
            }
        }
    }
}


int main() {
    APPS {
        APP(hw10a);
        APP(hw10b);
    }
    return 0;
}





