// #define COW_PATCH_FRAMERATE
// #define COW_PATCH_FRAMERATE_SLEEP
#include "include.cpp"

////////////////////////////////////////////////////////////////////////////////
// mini world //////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

struct OrbitCamera {
    double distance;
    double theta;
    double phi;
    double _angle_of_view;
};

mat4 orbit_camera_get_C(OrbitCamera *orbit) {
    return {}; // TODO return C
}

void orbit_camera_move(OrbitCamera *orbit) {
    if ((globals._mouse_owner == COW_MOUSE_OWNER_NONE) && globals.mouse_left_held) {
        // TODO overwrite orbit->theta
        // TODO overwrite orbit->phi
        // NOTE make sure you clamp phi or it will be possible to "pass over the north or south poles"

    }

    // TODO overwrite orbit->distance

}

// --------------------------------------------------

struct FPSCamera {
    vec3 origin;
    double theta;
    double phi;
    double _angle_of_view;
};

mat4 fps_camera_get_C(FPSCamera *human) {
    return {}; // TODO return C
}

void fps_camera_move(FPSCamera *human) {
    // // TODO overwrite human->origin
    // TODO hold W key to walking forward
    // TODO hold S key to walk backward
    // TODO hold A, D key to strafe left and right


    if (window_is_pointer_locked()) {
        // TODO overwrite human->theta
        // TODO overwrite human->phi
        // NOTE make sure you clamp phi or it will be possible to "look backward through your legs"

    }
}

// --------------------------------------------------

struct TrackingCamera {
    vec3 origin;
    vec3 *target;
    double _angle_of_view;
};

mat4 tracking_camera_get_C(TrackingCamera *track) {
    return {}; // TODO
}

// --------------------------------------------------

struct ArbitraryCamera {
    // note: using a mat4's to represent rotation is maybe not the best
    //       but it sure is convenient                                 
    vec3 origin;
    mat4 R;
    double _angle_of_view;
};

mat4 arbitrary_camera_get_C(ArbitraryCamera *plane) {
    return {}; // TODO
}

void arbitrary_camera_move(ArbitraryCamera *plane) {
    // // TODO overwrite plane->R
    // TODO hold W key to pitch down
    // TODO hold S key to pitch up
    // TODO hold A to yaw left
    // TODO hold D to yaw right
    // TODO hold J to roll left
    // TODO hold L to roll right


    // begin don't change this code
    plane->origin += transformVector(plane->R, V3(0, 0, -1));
    // end   don't change this code
}

// --------------------------------------------------

void draw_box_with_fake_shadows(mat4 PV, mat4 M, vec3 color) {
    library.soups.box.draw(PV * M, color, 3);
    { // shadows
        Soup3D shadow_box = library.soups.box;
        shadow_box.primitive = SOUP_QUADS;

        // TODO set this matrix to go from world coordinates of the object to world coordinates of the shadow
        // HINT my solution replaces four of the 0's with 1's (but which ones?! that is the question!)       
        // NOTE qthe (1, 1)-entry should be 0
        mat4 M_FakeShadow = {
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
        };

        vec4 shadow_color = V4(0, 0, 0, .2);

        // begin don't change this code
        {
            glEnable(GL_CULL_FACE);
            glCullFace(GL_BACK);
            shadow_box.draw(PV * M_FakeShadow * M, shadow_color);
            glDisable(GL_CULL_FACE);
        }
        // end   don't change this code
    }
}

void draw_world(mat4 P, mat4 V) {
    mat4 PV = P * V;

    { // little patch of grass
        double r = 100;
        vec3 verts[] = { { -r, 0, -r }, { -r, 0,  r }, {  r, 0,  r }, {  r, 0, -r } };
        soup_draw(PV, SOUP_QUADS, 4, verts, NULL, .7 * monokai.green);
    }
    { // trees
        draw_box_with_fake_shadows(PV, M4_Translation(40, 10, 25) * M4_Scaling(3, 10, 3), monokai.brown);
        draw_box_with_fake_shadows(PV, M4_Translation(40, 20, 25) * M4_Scaling(10, 10, 10), AVG(monokai.green, monokai.yellow));
        draw_box_with_fake_shadows(PV, M4_Translation(60, 15, 50) * M4_Scaling(3, 15, 3) * M4_RotationAboutYAxis(RAD(30)), monokai.brown);
        draw_box_with_fake_shadows(PV, M4_Translation(60, 30, 50) * M4_Scaling(10, 10, 10) * M4_RotationAboutYAxis(RAD(30)), monokai.green);
    }
    { // outer box
        vec3 vertex_colors[] = { // FORNOW
            .5 * monokai.black,  .5 * monokai.black,  .5 * monokai.blue,   .5 * monokai.purple,
            .5 * monokai.orange, .5 * monokai.yellow, .5 * monokai.black,  .5 * monokai.black,
            .5 * monokai.black,  .5 * monokai.black,  .5 * monokai.black,  .5 * monokai.black,
            .5 * monokai.brown,  .5 * monokai.gray,   .5 * monokai.purple, .5 * monokai.orange,
            .5 * monokai.black,  .5 * monokai.purple, .5 * monokai.orange, .5 * monokai.black,
            .5 * monokai.black,  .5 * monokai.yellow, .5 * monokai.blue,   .5 * monokai.black,
        };
        Soup3D outer_box = library.soups.box;
        outer_box.vertex_colors = vertex_colors;
        outer_box.primitive = SOUP_QUADS;
        outer_box.draw(PV * M4_Translation(0, 999.0, 0) * M4_Scaling(1000));
    }
}

// --------------------------------------------------

#define ORBIT_CAMERA 0
#define HUMAN_CAMERA 1
#define TRACK_CAMERA 2
#define PLANE_CAMERA 3
#define NO_ACTIVE_CAMERA -1
#define NUM_CAMERAS 4
#define INDEX2GRID(i) V2((i) % 2, (i) / 2) // 0 -> (0, 0); 1 -> (1, 0); ...

void hw04a() {

    OrbitCamera orbit = {};
    FPSCamera human = {};
    TrackingCamera track = {};
    ArbitraryCamera plane = {};
    int selected_camera = NO_ACTIVE_CAMERA;
    bool _initialized = false;
    bool draw_camera_cubes = true;
    window_set_clear_color(0.1, 0.1, 0.1, 1.0);

    while (cow_begin_frame()) {

        { // gui
            { // selected_camera
                int tmp = selected_camera;
                gui_slider("selected_camera", &selected_camera, -1, 3, 0, COW_KEY_TAB, true);

                // fornow
                if ((tmp != selected_camera)) {
                    window_pointer_unlock();
                }
            }
            gui_checkbox("draw_camera_cubes", &draw_camera_cubes, 'z');
        }

        { // init
            if (!_initialized || gui_button("reset", 'r')) {
                _initialized = true;
                orbit = { 200, RAD(15), RAD(-30), RAD(75) };
                human = { V3(0, 10, -20), RAD(180), 0, RAD(60), };
                track = { V3(-50, 50, 50), &plane.origin, RAD(45) };
                plane = { V3(0, 100, -500), M4_RotationAboutYAxis(RAD(180)), RAD(45) };
                selected_camera = NO_ACTIVE_CAMERA;
                window_pointer_unlock();
            }
        }

        { // update
            if (selected_camera == ORBIT_CAMERA) {
                orbit_camera_move(&orbit);
            } else if (selected_camera == PLANE_CAMERA) {
                arbitrary_camera_move(&plane);
            } else if (selected_camera == HUMAN_CAMERA) {
                fps_camera_move(&human);
            }
        }

        { // draw

            // matrices for all cameras
            mat4 C_array[] = {
                orbit_camera_get_C(&orbit),
                fps_camera_get_C(&human),
                tracking_camera_get_C(&track),
                arbitrary_camera_get_C(&plane), };
            mat4 P_array[] = {
                _window_get_P_perspective(orbit._angle_of_view),
                _window_get_P_perspective(human._angle_of_view),
                _window_get_P_perspective(track._angle_of_view),
                _window_get_P_perspective(plane._angle_of_view), };
            mat4 V_array[NUM_CAMERAS]; {
                for (int i = 0; i < NUM_CAMERAS; ++i) {
                    V_array[i] = inverse(C_array[i]);
                }
            }

            { // draw
                { // split screen bars
                    eso_begin(globals.Identity, SOUP_LINES, 6.0, false, true);
                    eso_color(0.0, 0.0, 0.0);
                    eso_vertex(-1.0, 0.0);
                    eso_vertex( 1.0, 0.0);
                    eso_vertex(0.0, -1.0);
                    eso_vertex(0.0,  1.0);
                    eso_end();
                }

                vec2 _window_radius = window_get_size() / 2;

                for (int camera_index = 0; camera_index < NUM_CAMERAS; ++camera_index) {
                    vec2 grid = INDEX2GRID(camera_index);
                    vec2 sgn = 2 * grid - V2(1.0);
                    text_draw(globals.Identity,
                            (char *) (
                                (camera_index == ORBIT_CAMERA) ? "orbit" :
                                (camera_index == HUMAN_CAMERA) ? "human" :
                                (camera_index == TRACK_CAMERA) ? "track" :
                                "plane"),
                            V2(-.96, -.9) + grid,
                            color_kelly(camera_index), 24, {}, true);

                    mat4 P = M4_Translation(.5 * sgn.x, .5 * sgn.y) * M4_Scaling(.5) * P_array[camera_index];
                    mat4 V = V_array[camera_index];
                    mat4 PV = P * V;

                    glEnable(GL_SCISSOR_TEST); // https://registry.khronos.org/OpenGL-Refpages/es2.0/xhtml/glScissor.xml
                    glScissor(int(grid[0] * _window_radius.x), int(grid[1] * _window_radius.y), int(_window_radius.x), int(_window_radius.y));
                    {
                        draw_world(P, V);

                        if (draw_camera_cubes) { // camera boxes and axes
                            for (int i = 0; i < NUM_CAMERAS; ++i) {
                                if (camera_index != i) { // don't draw ourself if "we are the camera"
                                    // note: this is very subtle; we're using C *as* M, because we're *drawing C itself*
                                    draw_box_with_fake_shadows(PV, C_array[i] * M4_Scaling(3), color_kelly(i));
                                    library.soups.axes.draw(PV * C_array[i] * M4_Scaling(10.0), color_kelly(i), 5.0);
                                    { // (potentially) extended -z axis
                                        eso_begin(PV, SOUP_LINES, 5.0);
                                        vec3 o = transformPoint(C_array[i], V3(0., 0., 0.));
                                        vec3 dir = normalized(transformVector(C_array[i], V3(0, 0, -1)));
                                        double L = (i == ORBIT_CAMERA) ? orbit.distance : (i == TRACK_CAMERA) ? norm(track.origin - *track.target) : 100;
                                        vec3 color = color_kelly(i);
                                        eso_color(color);
                                        eso_vertex(o);
                                        eso_color(color, .5);
                                        eso_vertex(o + L * dir);
                                        eso_end();
                                    }
                                }
                            }
                        }
                    }
                    glDisable(GL_SCISSOR_TEST);
                }
            }
        }

        { // gui
            { // picking
                if (globals._mouse_owner == COW_MOUSE_OWNER_NONE) {
                    int camera_hot = NO_ACTIVE_CAMERA;

                    auto gui_picker = [&](int camera_index) {
                        vec2 sgn = 2 * INDEX2GRID(camera_index) - V2(1.0);
                        if ((sgn.x * globals.mouse_position_NDC.x > 0) && (sgn.y * globals.mouse_position_NDC.y > 0)) {
                            camera_hot = camera_index;
                        }
                        if (camera_index == camera_hot) {
                            if (camera_index != selected_camera) {
                                eso_begin(globals.Identity, SOUP_QUADS);
                                eso_color(1.0, 1.0, 1.0, 0.3);
                                eso_vertex(0.0, 0.0);
                                eso_vertex(0.0, sgn.y);
                                eso_vertex(sgn.x, sgn.y);
                                eso_vertex(sgn.x, 0.0);
                                eso_end();
                            }
                        }
                        if (camera_index == selected_camera) {
                            eso_begin(globals.Identity, SOUP_LINE_LOOP, 24.0, false, true);
                            eso_color(color_kelly(camera_index));
                            eso_vertex(0.0, 0.0);
                            eso_vertex(0.0, sgn.y);
                            eso_vertex(sgn.x, sgn.y);
                            eso_vertex(sgn.x, 0.0);
                            eso_end();
                        }
                    };

                    gui_picker(HUMAN_CAMERA);
                    if ((selected_camera != HUMAN_CAMERA) || (!window_is_pointer_locked())) {
                        gui_picker(ORBIT_CAMERA);
                        gui_picker(TRACK_CAMERA);
                        gui_picker(PLANE_CAMERA);
                    }

                    if ((camera_hot != NO_ACTIVE_CAMERA) && globals.mouse_left_pressed) {
                        selected_camera = camera_hot;
                        if (selected_camera == HUMAN_CAMERA) {
                            window_pointer_lock();
                        }
                    }

                    if (globals.key_pressed[COW_KEY_ESCAPE]) {
                        if (selected_camera == HUMAN_CAMERA) { window_pointer_unlock(); }
                        selected_camera = NO_ACTIVE_CAMERA;
                    }
                }
            }

            { // tweaks
                if (selected_camera == ORBIT_CAMERA) {
                    gui_readout("distance", &orbit.distance);
                    gui_readout("theta", &orbit.theta);
                    gui_readout("phi", &orbit.phi);
                } else if (selected_camera == HUMAN_CAMERA) {
                    gui_readout("origin", &human.origin);
                    gui_readout("theta", &human.theta);
                    gui_readout("phi", &human.phi);
                    gui_printf("---");
                    gui_printf("click to lock pointer");
                    gui_printf("press Escape to unlock pointer");
                } else if (selected_camera == TRACK_CAMERA) {
                    if (gui_button("switch target", COW_KEY_SPACE)) {
                        if (track.target == &plane.origin) {
                            track.target = &human.origin;
                        } else {
                            track.target = &plane.origin;
                        }
                    }
                } else if (selected_camera == PLANE_CAMERA) {
                    gui_readout("altitude", &plane.origin.y);
                }
            }
        }

    }

}

////////////////////////////////////////////////////////////////////////////////
// the arcball /////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


void hw04b_compute_z(vec3 *v) {
    double sn_v = (squaredNorm(*v));
    v->z = (sn_v < 1)  ? sqrt(1 - sn_v) : 0;
}

void hw04b() {
    IndexedTriangleMesh3D bunny = library.meshes.bunny;
    Camera3D camera = { 5, RAD(45) };

    mat4 R_prev = globals.Identity;
    mat4 R_curr = globals.Identity;
    vec3 p_click = {};

    while (cow_begin_frame()) {

        if (globals.mouse_left_pressed) {
            // TODO

        } else if (globals.mouse_left_held) {
            // TODO (see reading)
            // NOTE i recommend _NOT_ calling normalized(u), since this will crash if |u| is 0
            // - instead, only update R_curr if (!IS_ZERO(squaredNorm(u)))                    

        } else if (globals.mouse_left_released) {
            // TODO

        }

        mat4 R_arcball = R_curr * R_prev; // don't change this line (unless, e.g., to debug something)

        // begin don't change this code
        {
            { // NOTE I am disabling built-in camera rotation and 2D offsets
                camera_move(&camera);
                camera.theta = 0;
                camera.phi = 0;
                camera.o_x = 0;
                camera.o_y = 0;
            }
            mat4 P = camera_get_P(&camera);
            mat4 V = camera_get_V(&camera);
            bunny.draw(P, V, R_arcball, monokai.blue);
            { // bespoke widget to visualize the arcball
                eso_begin(globals.Identity, SOUP_LINE_LOOP); { // ! we are supplying vertices _in NDC_ (how cool is that?)
                    eso_color(monokai.orange);
                    for (double theta = 0; theta < TAU - TINY_VAL; theta += (TAU) / 64) {
                        eso_vertex(e_theta(theta));
                    }
                } eso_end();
            }
        }
        // end   don't change this code

    }
}

int main() {
    APPS {
        APP(hw04a);
        APP(hw04b);
    }
    return 0;
}

