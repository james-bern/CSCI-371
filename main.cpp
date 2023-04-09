// #define COW_PATCH_FRAMERATE
// #define COW_PATCH_FRAMERATE_SLEEP
#include "include.cpp"

// NOTE you can drag the mass around with the mouse
//      (maybe pause -- press 'p' -- first) :)     

double h = .02;   // simulation timestep               
double L = 1;     // pendulum length                   
double g = -9.81; // gravitational acceleration        
double m = 1;     // mass (of point at end of pendulum)

vec2 hw9a_get_p(double theta) {
    // position of the pendulum's mass (end point)
    // convention: theta = 0 <-> straight down
    return L * V2(sin(theta), -cos(theta));
}

double hw9a_get_alpha(double theta) {
    // TODO angular acceleration
    return g / L * sin(theta);
};

double hw9a_get_PE(double theta) {
    // TODO potential energy
    double y = cos(theta) / L - L;
    return m * g * y;
}
double hw9a_get_KE(double omega) {
    // TODO kinetic energy
    double v = L * omega;
    return m * v * v / 2;
}
double hw9a_get_E(double theta, double omega) {
    // total energy
    return hw9a_get_PE(theta) + hw9a_get_KE(omega);
};

void hw9a() {
    Camera2D camera = { 5.35, 0, -.75 };
    bool paused = false;

    double theta = RAD(90);
    double omega = 0;

    const int TRACE_LENGTH = 128;
    StretchyBuffer<vec2> position_trace = {};
    vec3 trace_colors[TRACE_LENGTH]; {
        for (int node_i = 0; node_i < TRACE_LENGTH; ++node_i) {
            trace_colors[node_i] = color_plasma(double(node_i) / double(TRACE_LENGTH - 1));
        }
    }

    Plot plot = { 128 };
    plot_add_trace(&plot, 0.0, 16.0, monokai.blue);
    plot_add_trace(&plot, 0.0, 16.0, monokai.red);
    plot_add_trace(&plot, 0.0, 16.0, monokai.purple);

    enum Mode {
        MODE_EXPLICIT,
        MODE_SEMI_IMPLICIT,
        MODE_IMPLICIT,
        _NUM_MODES
    };
    int mode = 0;
    char *modes[_NUM_MODES];
    modes[MODE_EXPLICIT]      = "explicit euler";
    modes[MODE_SEMI_IMPLICIT] = "semi implicit euler";
    modes[MODE_IMPLICIT]      = "implicit euler";

    struct {
        bool hide_plots;
    } tweaks = {};

    while (cow_begin_frame()) {
        camera_move(&camera);
        mat4 PV = camera_get_PV(&camera);

        { // gui
            gui_slider("h", &h, .001, .1);
            gui_slider("L", &L, .1, 5);
            gui_slider("m", &m, .1, 10);
            gui_slider("g", &g, -24, 24);
            gui_readout("theta", &theta);
            gui_readout("omega", &omega);
            { // mode
                gui_slider("mode", &mode, 0, _NUM_MODES - 1, 'j', 'k', true);
                gui_printf("mode: %s", modes[mode]);
            }
            gui_checkbox("paused", &paused, 'p');
            { // reset
                static double _theta_0 = theta;
                static double _omega_0 = omega;
                if (gui_button("reset", 'r')) {
                    theta = _theta_0;
                    omega = _omega_0;
                    sbuff_free(&position_trace);
                }
            }
            gui_checkbox("hide_plots", &tweaks.hide_plots, 'b');
        }

        if (!paused) {
            // TODO step simulation state forward in time
            // NOTE please use hw9a_get_alpha(theta) here
            if (mode == MODE_EXPLICIT) {
                // TODO explicit euler (be careful!)
                double alpha = hw9a_get_alpha(theta);
                theta += h * omega;
                omega += h * alpha;
            } else if (mode == MODE_SEMI_IMPLICIT) {
                // TODO semi-implicit euler (be careful!)
                omega += h * hw9a_get_alpha(theta);
                theta += h * omega;
            } else if (mode == MODE_IMPLICIT) {
                // TODO implicit euler (harder; we will cover in Thursday's class)
                real theta_next = theta;
                real omega_next = omega;
                int iters = 0;
                while (true) {
                    iters++;
                    real theta_tmp = theta + h * omega_next;
                    real omega_tmp = omega + h * hw9a_get_alpha(theta_next);
                    if ((pow(theta_tmp - theta_next, 2) + pow(omega_tmp - omega_next, 2)) < .001 || iters > 100) {
                        theta = theta_tmp;
                        omega = omega_tmp;
                        break;
                    }
                    theta_next = theta_tmp;
                    omega_next = omega_tmp;
                }
            }
        }

        vec2 p = hw9a_get_p(theta);

        { // drag
            if (widget_drag(PV, 1, &p)) {
                p = L * normalized(p);
                omega = 0;
                theta = RAD(90) + atan2(p.y, p.x);
            }
        }

        { // draw
            { // position_trace
                if (position_trace.length < TRACE_LENGTH) {
                    sbuff_push_back(&position_trace, p);
                } else {
                    memmove(position_trace.data, position_trace.data + 1, (TRACE_LENGTH - 1) * sizeof(vec2));
                    position_trace.data[TRACE_LENGTH - 1] = p;
                }
                soup_draw(PV, SOUP_LINE_STRIP, position_trace.length, position_trace.data, trace_colors);
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
            }
            { // energy plot
                if (!tweaks.hide_plots) {
                    mat4 PV_plot = M4_Translation(.3, -.3) * M4_Scaling(.5, .5 * _window_get_aspect());
                    if (!paused) {
                        plot_data_point(&plot, 0, hw9a_get_PE(theta));
                        plot_data_point(&plot, 1, hw9a_get_KE(omega));
                        plot_data_point(&plot, 2, hw9a_get_E(theta, omega));
                    }
                    plot_draw(&plot, PV_plot);
                    text_draw(PV_plot, "time", { 1.05, 0.0 });
                    text_draw(PV_plot, "energy", { 0.05, 1.0 });
                }
            }
        }
    }
}




void hw9b() {
    // // skeleton                                                               
    // bone_lengths          -- the length of the bones                          
    // bone_relative_angles  -- the angle between a given bone and the next bone 
    // bone_rest_positions   -- the origin of each bone in the initial pose      
    // * NOTE You may assume that the rest positions all lie along the x-axis.   
    // skeleton_total_length -- the lenth of the entire skeleton in its rest pose
    const int NUM_BONES = 4;
    double bone_lengths[NUM_BONES];
    double bone_relative_angles[NUM_BONES] = {};
    vec2 bone_rest_positions[NUM_BONES + 1] = {};
    double skeleton_total_length = 0;
    {
        for (int bone_i = 0; bone_i < NUM_BONES; ++bone_i) {
            bone_lengths[bone_i] = LERP(double(bone_i) / MAX(1, NUM_BONES - 1), 2.5, 1);
            bone_rest_positions[bone_i + 1] = bone_rest_positions[bone_i] + V2(bone_lengths[bone_i], 0);
        }
        skeleton_total_length = bone_rest_positions[NUM_BONES].x;
    }

    enum Mode {
        MODE_RIGID,
        MODE_SMOOTH,
        _NUM_MODES,
    };
    char *modes[_NUM_MODES];
    modes[MODE_RIGID]  = "rigid";
    modes[MODE_SMOOTH] = "smooth";

    // // skin                                                                           
    // node_rest_positions -- the positions of all nodes in the rest pose                
    // weights             -- for a given skinning mode, the binding weights of each node
    const int NUM_NODES = 64; STATIC_ASSERT(NUM_NODES % 2 == 0);
    vec2 node_rest_positions[NUM_NODES] = {};
    double weights[_NUM_MODES][NUM_NODES][NUM_BONES] = {}; // weights[mode][node_i][bone_i]
    {
        {
            // rest pose as a loop
            // e.g., for 10 nodes:
            // 9 8 7 6 5          
            // 0 1 2 3 4          
            int k = 0;
            for (int sign = -1; sign <= 1; sign += 2) {
                for (int node_i = 0; node_i < NUM_NODES / 2; ++node_i) {
                    double f = double(node_i) / double(NUM_NODES / 2 - 1);
                    if (sign > 0) { f = 1 - f; }
                    node_rest_positions[k++] = V2(skeleton_total_length * f, .3 * sign);
                }
            }
        }

        // weights
        for (int node_i = 0; node_i < NUM_NODES / 2; ++node_i) {
            for (int mode = 0; mode < _NUM_MODES; ++mode) {
                // calculate the node_i-th node's weights (not normlized)
                double w[NUM_BONES] = {};
                {

                    if (mode == MODE_RIGID) {
                        // idea: bind the node entirely to the closest bone
                        // NOTE I did this for you                         
                        double X_i = node_rest_positions[node_i].x;
                        for (int bone_i = 0; bone_i < NUM_BONES; ++bone_i) {
                            if (X_i < bone_rest_positions[bone_i + 1].x + TINY_VAL) {
                                w[bone_i] = 1;
                                break;
                            }
                        }
                    }

                    else if (mode == MODE_SMOOTH) {
                        // idea: smoothly blend the weights between nearby bones        
                        // (many possible ways to implement this)                       
                        // TODO (harder; we will cover one approach in Thursday's class)
                        double X_i = node_rest_positions[node_i].x;
                        double c_kernel = X_i;
                        double r_kernel = .6;
                        double L_kernel = MAX(0, c_kernel - r_kernel);
                        double R_kernel = MIN(c_kernel + r_kernel, skeleton_total_length);
                        for (int bone_i = 0; bone_i < NUM_BONES; ++bone_i) {
                            double L_bone = bone_rest_positions[bone_i].x;
                            double R_bone = bone_rest_positions[bone_i + 1].x;
                            if (L_bone - TINY_VAL < L_kernel && R_bone + TINY_VAL > R_kernel) {
                                w[bone_i] += R_kernel - L_kernel;
                            } else if (IS_BETWEEN(L_bone, L_kernel, R_kernel)) {
                                w[bone_i] += R_kernel - L_bone;
                            } else if (IS_BETWEEN(R_bone, L_kernel, R_kernel)) {
                                w[bone_i] += R_bone - L_kernel;
                            }
                        }
                    }
                }

                { // normalize
                    double sum_w = 0; {
                        for (int bone_i = 0; bone_i < NUM_BONES; ++bone_i) {
                            sum_w += w[bone_i];
                        }
                        for (int bone_i = 0; bone_i < NUM_BONES; ++bone_i) {
                            w[bone_i] /= sum_w;
                        }
                    }
                }

                { // write w into 3D array
                    for (int bone_i = 0; bone_i < NUM_BONES; ++bone_i) {
                        weights[mode][node_i][bone_i] = weights[mode][NUM_NODES - 1 - node_i][bone_i] = w[bone_i];
                    }
                }
            }
        }
    }

    struct {
        bool draw_character;
        bool hide_skeleton;
        bool hide_nodes;
        bool draw_rest_pose;
        bool hide_plots;
    } tweaks = {};


    // debug play
    bool debug_play = false;
    double debug_time = 0;

    // keyframing system
    const int MAX_KEYFRAMES = 1024;
    double keyframes[MAX_KEYFRAMES][NUM_BONES] = {};
    int num_keyframes = 0;
    bool tween_keyframes = false;
    const int FRAMES_PER_KEYFRAME = 32;
    int frame = 0;

    Camera2D camera = { 10, 0, 2 };
    int mode = 0;

    while (cow_begin_frame()) {
        camera_move(&camera);
        mat4 PV = camera_get_PV(&camera);

        { // gui
            for (int bone_i = 0; bone_i < NUM_BONES; ++bone_i) {
                char buffer[24];
                sprintf(buffer, "theta_%d", bone_i);
                gui_slider(buffer, bone_relative_angles + bone_i, -PI, PI, true);
            }
            gui_slider("mode", &mode, 0, _NUM_MODES - 1, 'j', 'k', true);
            gui_printf("mode: %s", modes[mode]);
            if (!tweaks.hide_plots) {
                mat4 PV_plot = M4_Translation(0, .25) * M4_Scaling(.59, .5 * .5 * _window_get_aspect());
                static vec2 trace[NUM_NODES / 2];
                {
                    vec2 axes[] = { { 1, 0 }, { 0, 0 }, { 0, 1 } };
                    soup_draw(PV_plot, SOUP_LINE_STRIP, 3, axes, NULL, monokai.gray);
                    text_draw(PV_plot, "x_rest", axes[0]);
                    text_draw(PV_plot, "weight", axes[2], V3(1, 1, 1), 0, { 0, -24 });
                }
                for (int bone_i = 0; bone_i < NUM_BONES; ++bone_i) {
                    for (int node_i = 0; node_i < NUM_NODES / 2; ++node_i) {
                        trace[node_i] = { real(node_i) / real(NUM_NODES / 2 - 1), weights[mode][node_i][bone_i] };
                    }
                    soup_draw(PV_plot, SOUP_LINE_STRIP, NUM_NODES / 2, trace, NULL, color_kelly(bone_i));
                }
            }
            gui_checkbox("debug_play", &debug_play, 'p');
            if (gui_button("save keyframe", 's') && num_keyframes < MAX_KEYFRAMES) {
                for (int bone_i = 0; bone_i < NUM_BONES; ++bone_i) {
                    keyframes[num_keyframes][bone_i] = bone_relative_angles[bone_i];
                }
                ++num_keyframes;
            }
            gui_checkbox("tween keyframes", &tween_keyframes, 't');
            gui_readout("num_keyframes", &num_keyframes);
            { // reset
                static double _bone_relative_angles_0[NUM_BONES];
                do_once { memcpy(_bone_relative_angles_0, bone_relative_angles, sizeof(bone_relative_angles)); };
                if (gui_button("reset", 'r')) {
                    memcpy(bone_relative_angles, _bone_relative_angles_0, sizeof(bone_relative_angles));

                    debug_time = 0;

                    memset(keyframes, 0, sizeof(keyframes));
                    num_keyframes = 0;
                    tween_keyframes = false;
                    frame = 0;
                }
            }
            gui_checkbox("draw_character", &tweaks.draw_character, 'z');
            gui_checkbox("draw_rest_pose", &tweaks.draw_rest_pose, 'v');
            gui_checkbox("hide_skeleton", &tweaks.hide_skeleton, 'x');
            gui_checkbox("hide_nodes", &tweaks.hide_nodes, 'c');
            gui_checkbox("hide_plots", &tweaks.hide_plots, 'b');
        }

        // aniamte the skeleton (i.e., set bone_relative_angles using whatever method you like)
        // METHOD 1: press 'p' to play a sinusoidal trajectory                                 
        // METHOD 2: press 's' to save keyframes; press 't' to play them back with lerp        
        {
            if (tween_keyframes && num_keyframes > 0) {
                int A = (frame / FRAMES_PER_KEYFRAME) % num_keyframes;
                int B = (A + 1) % num_keyframes;
                double t = double(frame % FRAMES_PER_KEYFRAME) / FRAMES_PER_KEYFRAME;
                for (int bone_i = 0; bone_i < NUM_BONES; ++bone_i) {
                    bone_relative_angles[bone_i] = LERP(
                            t,
                            keyframes[A][bone_i],
                            keyframes[B][bone_i]);
                }
                ++frame;
            } else if (debug_play) {
                for (int bone_i = 1; bone_i < NUM_BONES; ++bone_i) {
                    bone_relative_angles[bone_i] += .02 * cos((bone_i + 1) * debug_time / 2);
                }
                debug_time += .0167;
            }
        }

        // bone_absolute_angles -- the angle of each bone in world coordinates
        // NOTE feel free to use these for forward kinematics and skinning!   
        double bone_absolute_angles[NUM_BONES] = {};
        {
            for (int bone_i = 0; bone_i < NUM_BONES; ++bone_i) {
                if (bone_i > 0) {
                    bone_absolute_angles[bone_i] = bone_absolute_angles[bone_i - 1];
                }
                bone_absolute_angles[bone_i] += bone_relative_angles[bone_i];
            }
        }

        // // forward kinematics                                          
        // TODO compute bone_current_positions (the skeleton will show up)
        vec2 bone_current_positions[NUM_BONES + 1] = {};
        {
            for (int bone_i = 0; bone_i < NUM_BONES + 1; ++bone_i) {
                if (bone_i > 0) {
                    bone_current_positions[bone_i] = bone_current_positions[bone_i - 1] \
                                                     + bone_lengths[bone_i - 1] * e_theta(bone_absolute_angles[bone_i - 1]);
                }
            }
        }


        // // skinning                                             
        // TODO compute node_curr_positions (the skin will show up)
        vec2 node_curr_positions[NUM_NODES] = {};
        {
            for (int node_i = 0; node_i < NUM_NODES; ++node_i) {
                vec2 P = node_rest_positions[node_i];
                for (int bone_i = 0; bone_i < NUM_BONES; ++bone_i) {
                    vec2 O = bone_rest_positions[bone_i];
                    vec2 S = P - O;

                    double phi = bone_absolute_angles[bone_i];
                    vec2 p = bone_current_positions[bone_i];
                    node_curr_positions[node_i] += weights[mode][node_i][bone_i] * (rotated(S, phi) + p);
                }
            }
        }

        { // draw
            if (!tweaks.hide_skeleton) { // draw
                eso_begin(PV, GL_LINES); {
                    for (int bone_i = 0; bone_i < NUM_BONES; ++bone_i) {
                        eso_color(color_kelly(bone_i));
                        eso_vertex(bone_current_positions[bone_i]);
                        eso_vertex(bone_current_positions[bone_i + 1]);
                    }
                } eso_end();
                soup_draw(PV, SOUP_POINTS, NUM_BONES + 1, bone_current_positions, NULL, monokai.white);
            }

            if (tweaks.draw_rest_pose) {
                soup_draw(PV, SOUP_POINTS, NUM_BONES + 1, bone_rest_positions, NULL, monokai.blue);
                soup_draw(PV, SOUP_LINE_STRIP, NUM_BONES + 1, bone_rest_positions, NULL, monokai.blue);
                soup_draw(PV, SOUP_LINE_LOOP, NUM_NODES, node_rest_positions, NULL, .5 * monokai.blue);
                soup_draw(PV, SOUP_POINTS, NUM_NODES, node_rest_positions, NULL, monokai.blue, 8.0);
            }

            if (!tweaks.hide_nodes) {
                soup_draw(PV, SOUP_LINE_LOOP, NUM_NODES, node_curr_positions, NULL, monokai.gray);
                soup_draw(PV, SOUP_POINTS, NUM_NODES, node_curr_positions, NULL, monokai.white, 8.0);
            }

            if (tweaks.draw_character) {
                const int NUM_QUADS = (NUM_NODES / 2 - 1);
                static vec2 vertex_positions[4 * NUM_QUADS];
                static vec3 vertex_colors[4 * NUM_QUADS];
                for (int node_i = 0; node_i < NUM_QUADS; ++node_i) {
                    vertex_positions[4 * node_i + 0] = node_curr_positions[node_i];
                    vertex_positions[4 * node_i + 1] = node_curr_positions[node_i + 1];
                    vertex_positions[4 * node_i + 2] = node_curr_positions[NUM_NODES - 1 - (node_i + 1)];
                    vertex_positions[4 * node_i + 3] = node_curr_positions[NUM_NODES - 1 - (node_i)];
                    double f = double(node_i) / (NUM_QUADS - 1);
                    vec3 color = color_rainbow_swirl(f);
                    vertex_colors[4 * node_i + 0] = color;
                    vertex_colors[4 * node_i + 1] = color;
                    vertex_colors[4 * node_i + 2] = color;
                    vertex_colors[4 * node_i + 3] = color;
                }
                soup_draw(PV, SOUP_QUADS, 4 * NUM_QUADS, vertex_positions, vertex_colors);
            }
        }
    }
}


int main() {
    APPS {
        APP(hw9a);
        APP(hw9b);
    }
    return 0;
}



