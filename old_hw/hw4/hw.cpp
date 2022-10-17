#define COW_CRASH_ON_FLOATING_POINT_EXCEPTIONS
#define _CRT_SECURE_NO_WARNINGS
#include "snail.cpp"
#include "cow.cpp"
#include "_cow_supplement.cpp" // my sol'n to load_basic_mesh(...) and transform_vertex_positions_to_double_unit_box(...)

////////////////////////////////////////////////////////////////////////////////
// hw4a() //////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// NOTE I recommend spending 10-20 minutes playing with this app after you get it working.

void hw4a() {
    init();

    BasicTriangleMesh3D basic_box = load_basic_mesh("data_basic_box", true);
    BasicTriangleMesh3D basic_bunny = load_basic_mesh("data_basic_bunny", true);
    BasicTriangleMesh3D basic_teapot = load_basic_mesh("data_basic_teapot", true);
    int mesh = 0;

    double height_of_camera_field_of_view_at_origin = 8; // H in figure below

    Camera3D camera = {};
    camera.distance_to_origin = 15; // D in figure below
    camera.angle_of_view = RAD(45); // theta in figure below

    while (begin_frame()) {

        imgui_slider("height_of_camera_field_of_view_at_origin", &height_of_camera_field_of_view_at_origin, 1, 20);
        { // angle_of_view slider in degrees
            int tmp = int(round(DEG(camera.angle_of_view)));
            imgui_slider("angle_of_view (degrees)", &tmp, 1, 179);
            camera.angle_of_view = RAD(tmp);
        }
        { // NOTE I am disabling all built-in camera controls besides theta, phi
            input._mouse_wheel_offset = 0;
            camera_move(&camera);
            camera._o_x = 0;
            camera._o_y = 0;
        }
        imgui_slider("mesh", &mesh, 0, 2, 'j', 'k', true);
        imgui_readout("camera", &camera);

        { // TODO set camera.distance to origin

            // // orbit camera setup               
            //                                     
            //                      +   -  -  -   ^
            //                  -   |             |
            //               -      |             |
            //            -         |             |
            //         -            |             |
            // camera-) theta       origin        H
            //         -            |             |
            //            -         |             |
            //               -      |             |
            //                  -   |             |
            //                      +   -  -  -   V
            //                                     
            //       <------D------->              

            // HINT Your answer should involve dividing by 2 twice.

            // camera.distance_to_origin = ...; // what is D as a function of theta and H?
        }

        mat4 PV = camera_get_PV(&camera);
        if (mesh == 0) {
            basic_draw(TRIANGLE_MESH, PV, basic_box, monokai.green, 3, false, LERP(.4, monokai.green, monokai.white));
        } else if (mesh == 1) {
            basic_draw(TRIANGLE_MESH, PV, basic_bunny, monokai.blue, 3, false, LERP(.4, monokai.blue, monokai.white));
        } else {
            basic_draw(TRIANGLE_MESH, PV, basic_teapot, monokai.red, 3, false, LERP(.4, monokai.red, monokai.white));
        }

    }
}

////////////////////////////////////////////////////////////////////////////////
// hw4b() //////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// // you may want these functions                      
//                                                      
// vec3 transformPoint(mat4 M, vec3 p);                 
// vec3 transformVector(mat4 M, vec3 v);                
//                                                      
// mat4 RotationX(double angle_in_radians);             
// mat4 RotationY(double angle_in_radians);             
// mat4 RotationZ(double angle_in_radians);             
// mat4 Translation(vec3 t);                            
// mat4 Translation(double t_x, double t_y, double t_z);
// mat4 Scaling(vec3 s);                                
// mat4 Scaling(double s_x, double s_y, double s_z);    
//                                                      
// mat4 xyzo2mat4(vec3 x, vec3 y, vec3 z, vec3 o);      
//                                                      
// bool IS_ZERO(double a); // whether a is approx. zero 
// double norm(vecX v); // length of v                  
#if 0
vec3 ds = transformVector(RotationY(human->theta), V3(0, 0, -1));
#endif

// // you may want this input state                                                                              
//                                                                                                               
// input.key_held['w']       // (bool) whether W key is currently held                                           
// input.mouse_left_held     // (bool) whether mouse left button is currently held                               
// input._mouse_wheel_offset // (double) how far mouse wheel moved since last frame                              
// input._mouse_dx_NDC       // (double) how far mouse cursor moved in x since last frame mapped to range [-1, 1]
// input._mouse_dy_NDC       // (double) how far mouse cursor moved in y since last frame mapped to range [-1, 1]
#if 0
if (input.key_held['w']) {
    human->origin += ds;
}
#endif

// --------------------------------------------------

// feel free to add your own tweaks :)

struct {
    bool draw_axes = true;
    bool kelly_color_axes = false;
    bool extend_negative_z_axis = true;
    bool human_AI = false;
    double axes_scale = 10.;
    bool draw_labels = true;
    bool draw_fake_shadows = true;
} tweaks;

// --------------------------------------------------

struct OrbitCamera {
    double distance;
    double angle_of_view;
    double theta;
    double phi;
};

mat4 orbit_camera_get_C(OrbitCamera *orbit) { FORNOW_UNUSED(orbit);
    // TODO (see slides)
    return Identity4x4;
}

void orbit_camera_move(OrbitCamera *orbit) { FORNOW_UNUSED(orbit);
    if (input.mouse_left_held) {
        // TODO overwrite orbit->theta
        // TODO overwrite orbit->phi (make sure you clamp phi or it will be possible to "pass over the north or south poles")
    }
    // TODO overwrite orbit->distance (when user moves the mouse wheel)
}

// --------------------------------------------------

struct FPSCamera {
    vec3 origin;
    double angle_of_view;
    double theta;
    double phi;
};

mat4 fps_camera_get_C(FPSCamera *human) { FORNOW_UNUSED(human);
    // TODO (see slides)
    return Identity4x4;
}

void fps_camera_move(FPSCamera *human) { FORNOW_UNUSED(human);
    // // TODO overwrite human->origin
    // TODO hold W key to walking forward
    // TODO hold S key to walk backward
    // TODO hold A, D key to strafe left and right

    if (glfwGetInputMode(window, GLFW_CURSOR) == GLFW_CURSOR_DISABLED) { // pointer lock
        // TODO overwrite human->theta
        // TODO overwrite human->phi (make sure you clamp phi or it will be possible to "look backward through your legs")
    }
}

// --------------------------------------------------

struct TrackingCamera {
    vec3 origin;
    double angle_of_view;
    vec3 *target;
};

mat4 tracking_camera_get_C(TrackingCamera *track) { FORNOW_UNUSED(track);
    // TODO (see slides)
    // HINT: mat4 xyzo2mat4(vec3 x, vec3 y, vec3 z, vec3 o);
    return Identity4x4;
}

// --------------------------------------------------

struct ArbitraryCamera {
    // note: using a mat4's to represent rotation is maybe not the best
    //       but it sure is convenient                                 
    vec3 origin;
    double angle_of_view;
    mat4 R;
};

mat4 arbitrary_camera_get_C(ArbitraryCamera *plane) { FORNOW_UNUSED(plane);
    return Identity4x4;
}

void arbitrary_camera_move(ArbitraryCamera *plane) {
    // // TODO overwrite plane->R
    // TODO hold W key to pitch down
    // TODO hold S key to pitch up
    // TODO hold A to yaw left
    // TODO hold D to yaw right
    // TODO hold J to roll left
    // TODO hold L to roll right

    // begin you are _not_ allowed to change this line or work around it
    plane->origin += 4 * transformVector(plane->R, V3(0, 0, -1));
    // end   you are _not_ allowed to change this line or work around it
}

// --------------------------------------------------

vec3 basic_box_vertex_positions[] = {
    {-1,-1,1}, {-1,-1,-1}, {-1,1,-1}, {-1,1,1}, // left
    {1,1,1},{1,1,-1},{1,-1,-1},{1,-1,1},        // right
    {1,-1,1},{1,-1,-1},{-1,-1,-1},{-1,-1,1},    // bottom
    {-1,1,1}, {-1,1,-1}, {1,1,-1}, {1,1,1},     // top
    {-1,1,-1}, {-1,-1,-1}, {1,-1,-1}, {1,1,-1}, // back
    {1,1,1},{1,-1,1},{-1,-1,1},{-1,1,1},        // front
};

void draw_basic_box_with_fake_shadows(mat4 PV, mat4 M, vec3 color) {
    // draw the box itself
    basic_draw(QUAD_MESH, PV * M, 24, basic_box_vertex_positions, color, 3, false, .5 * color);

    if (tweaks.draw_fake_shadows) {
        // TODO set this matrix to go from world coordinates of the object to world coordinates of the shadow
        // HINT my solution has 4 non-zero entries
        mat4 M_FakeShadow = {
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
        };

        // TODO set this vector V4(red, green, blue, alpha)
        vec4 shadow_color = V4(0, 0, 0, 0);

        // begin don't change this code
        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);
        basic_draw(QUADS, PV * M_FakeShadow * M, 24, basic_box_vertex_positions, shadow_color);
        glDisable(GL_CULL_FACE);
        // end   don't change this code
    }
}

// --------------------------------------------------

enum Mode {               Orbit ,  Human ,  Track ,  Plane, NUM_MODES };
char *MODE_STRINGS[] = { "Orbit", "Human", "Track", "Plane" };

void hw4b() {
    init();

    OrbitCamera orbit = { 200, RAD(75), RAD(15), RAD(-30) };
    FPSCamera human = { V3(0, 10, -20), RAD(60), RAD(180), 0 };
    TrackingCamera track = { V3(-50, 50, 50), RAD(45), &human.origin };
    ArbitraryCamera plane = { V3(0, 100, -500), RAD(45), RotationY(RAD(180)) };

    Mode mode = Orbit;

    bool gui = true; // SET THIS FALSE IF RUNNING SLOW

    while (begin_frame()) {

        imgui_checkbox("", &gui);
        if (gui) { // gui
            { // tweaks
                imgui_checkbox("human_AI", &tweaks.human_AI, 'z');
                imgui_checkbox("draw_axes", &tweaks.draw_axes, 'x');
                if (tweaks.draw_axes) {
                    imgui_checkbox("extend_negative_z_axis", &tweaks.extend_negative_z_axis, 'c');
                    imgui_checkbox("kelly_color_axes", &tweaks.kelly_color_axes, 'v');
                    imgui_slider("axes_scale", &tweaks.axes_scale, 1., 100.);
                }
                imgui_checkbox("draw_labels", &tweaks.draw_labels, 'b');
                imgui_checkbox("draw_fake_shadows", &tweaks.draw_fake_shadows, 'g');
            }

            { // mode, pointer lock
                int prev_mode = mode;
                imgui_slider("mode", (int *) &mode, 0, NUM_MODES - 1, 0, KEY_TAB, true);
                for (int i = 0; i < NUM_MODES; ++i) {
                    if (input.key_pressed['0' + i]) {
                        mode = (Mode) i;
                    }
                }
                _imgui_printf("mode = %s", MODE_STRINGS[mode]);
                if (mode != prev_mode) {
                    if (mode != Human) {
                        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
                    }
                }
                if (mode == Human) {
                    if (widget_active_widget_ID == 0 && input.mouse_left_pressed) {
                        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
                    } else if (input.key_pressed[KEY_ESCAPE]) {
                        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
                    }
                }
            }

            if (mode == Plane) { imgui_readout("y", &plane.origin.y); }

            { // angle_of_view slider in degrees
                double *angle_of_view = (mode == Orbit) ? &orbit.angle_of_view : (mode == Human) ? &human.angle_of_view : (mode == Track) ? &track.angle_of_view : &plane.angle_of_view;
                double tmp = DEG(*angle_of_view);
                imgui_slider("angle_of_view (degrees)", &tmp, 1, 179);
                *angle_of_view = RAD(tmp);
            }
        }

        { // movement
            if (mode == Orbit) {
                orbit_camera_move(&orbit);
            } else if (mode == Plane) {
                arbitrary_camera_move(&plane);
            } else if (mode == Human && !tweaks.human_AI) {
                fps_camera_move(&human);
            }
            if (tweaks.human_AI) { // walk to random position
                static bool started;
                static int elapsed_frames;
                static int turn_sign;
                static int turn_num_frames;
                static int walk_num_frames;
                double walk_speed = .5;
                double turn_speed = .03;
                if (!started) {
                    started = true;
                    elapsed_frames = 0;
                    vec2 target = V2(util_random_double(-100, 100), util_random_double(-100, 100));
                    vec2 v = target - V2(human.origin.x, human.origin.z);
                    double theta = atan2(v.x, v.y) + RAD(180); // ??
                    double dtheta = WRAP(theta - human.theta, -PI, PI);
                    turn_sign = SGN(dtheta);
                    double mag_dtheta = ABS(dtheta);
                    double mag_ds = norm(v);
                    turn_num_frames = int(mag_dtheta / turn_speed);
                    walk_num_frames = int(mag_ds / walk_speed);
                }
                if (elapsed_frames < turn_num_frames) {
                    human.theta += turn_sign * turn_speed;
                } else if (elapsed_frames < turn_num_frames + walk_num_frames) {
                    human.origin += walk_speed * transformVector(RotationY(human.theta), V3(0, 0, -1));
                } else {
                    started = false;
                }
                ++elapsed_frames;
            }
        }


        { // draw

            // all cameras
            mat4 C[] = { orbit_camera_get_C(&orbit), fps_camera_get_C(&human), tracking_camera_get_C(&track), arbitrary_camera_get_C(&plane), };
            mat4 P[] = { tform_get_P_perspective(orbit.angle_of_view), tform_get_P_perspective(human.angle_of_view), tform_get_P_perspective(track.angle_of_view), tform_get_P_perspective(plane.angle_of_view), };
            mat4 V[NUM_MODES]; {
                for (int i = 0; i < NUM_MODES; ++i) {
                    V[i] = inverse(C[i]);
                }
            }

            // the camera actually rendering the scene
            mat4 P_mode = P[mode];
            mat4 V_mode = V[mode];
            mat4 PV_mode = P_mode * V_mode;

            { // little patch of grass
                double r = 100;
                vec3 verts[] = { { -r, 0, -r }, { -r, 0,  r }, {  r, 0,  r }, {  r, 0, -r } };
                vec3 colors[] = { monokai.green, AVG(monokai.green, monokai.blue), AVG(monokai.green, monokai.white), AVG(monokai.green, monokai.red) };
                basic_draw(QUADS, PV_mode, 4, verts, colors);
            }
            { // trees
                draw_basic_box_with_fake_shadows(PV_mode, Translation(40, 10, 0) * Scaling(3, 10, 3), monokai.brown);
                draw_basic_box_with_fake_shadows(PV_mode, Translation(40, 20, 0) * Scaling(10, 10, 10), AVG(monokai.green, monokai.yellow));
                draw_basic_box_with_fake_shadows(PV_mode, Translation(60, 15, -25) * Scaling(3, 15, 3) * RotationY(RAD(30)), monokai.brown);
                draw_basic_box_with_fake_shadows(PV_mode, Translation(60, 30, -25) * Scaling(10, 10, 10) * RotationY(RAD(30)), monokai.green);
            }
            { // sky (and ground) box
                vec3 vertex_colors[] = { // fornow: had to match orderings of the walls by trial and error
                    .5 * monokai.blue, .5 * monokai.blue, .1 * monokai.blue, .1 * monokai.blue, // left
                    .1 * monokai.blue, .1 * monokai.blue, .5 * monokai.blue, .5 * monokai.blue, // right
                    monokai.brown, monokai.gray, monokai.purple, monokai.orange,                // ground
                    .1 * monokai.blue, .1 * monokai.blue, .1 * monokai.blue, .1 * monokai.blue, // sky
                    .1 * monokai.blue, .5 * monokai.blue, .5 * monokai.blue, .1 * monokai.blue, // back
                    .1 * monokai.blue, .5 * monokai.blue, .5 * monokai.blue, .1 * monokai.blue, // front
                };
                basic_draw(QUADS, PV_mode * Translation(0, 999.0, 0) * Scaling(1000), 24, basic_box_vertex_positions, vertex_colors);
            }
            { // objects, axes, labels
                BasicMesh axes; {
                    BasicMesh _basic_axes; {
                        _basic_axes = meshlib.basic_axes;
                        _basic_axes.vertex_colors = NULL;
                    }
                    axes = (tweaks.kelly_color_axes) ? _basic_axes : meshlib.basic_axes;
                }
                for (int i = 0; i < NUM_MODES; ++i) {
                    if (mode != i) { // down draw ourself if "we are the camera"
                        draw_basic_box_with_fake_shadows(PV_mode, C[i] * Scaling(3), color_get_kelly(i));
                        if (tweaks.draw_labels) {
                            basic_text(PV_mode, MODE_STRINGS[i], transformPoint(C[i], V3(0., 0., 0.)), monokai.black);
                        }
                        if (tweaks.draw_axes) {
                            basic_draw(PV_mode * C[i] * Scaling(tweaks.axes_scale), axes, color_get_kelly(i));
                            if (tweaks.extend_negative_z_axis) {
                                gl_PV(PV_mode);
                                gl_begin(LINES);
                                vec3 o = transformPoint(C[i], V3(0., 0., 0.));
                                vec3 dir = normalized(transformVector(C[i], V3(0, 0, -1)));
                                double L = (i == Orbit) ? orbit.distance : (i == Track) ? norm(track.origin - *track.target) : 100;
                                vec3 color = color_get_kelly(i);
                                gl_color(color);
                                gl_vertex(o);
                                gl_color(color, .5);
                                gl_vertex(o + L * dir);
                                gl_end();
                            }
                        }
                    }
                }
            }
            { // TODO boundy castle, clouds
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
// hw4c() //////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// // you may want these functions        
// mat4 Rotation(vec3 axis, double angle);
// double squaredNorm(vec3 v); // |v|^2   
// double norm(vec3 v); // |v|            

// NOTE Identity4x4 is a 4x4 identity matrix you can copy
#if 0
R_curr = Identity4x4;
#endif

// // you may want this input state                                                                        
//                                                                                                         
// input.mouse_left_pressed  // (bool) whether mouse left button was just pressed                          
// input.mouse_left_held     // (bool) whether mouse left button is currently held                         
// input.mouse_left_released // (bool) whether mouse left button was just released                         
// input._mouse_x_NDC  // (double) mouse cursor x position mapped to range [-1, 1]                         
// input._mouse_y_NDC  // (double) mouse cursor y position mapped to range [-1, 1]                         

void hw4c() {
    init();

    FancyTriangleMesh3D fancy_bunny = load_fancy_mesh("data_fancy_bunny", true, true);
    BasicTriangleMesh3D basic_box = load_basic_mesh("data_basic_box", true);
    int mesh = 0;
    Camera3D camera = { 5, RAD(45) };

    mat4 R_prev = Identity4x4;
    mat4 R_curr = Identity4x4;
    vec3 p_click = {}; FORNOW_UNUSED(p_click);

    while (begin_frame()) {

        if (input.mouse_left_pressed) {
            // TODO (see reading)
        } else if (input.mouse_left_held) {
            // TODO (see reading)

            // NOTE: I recommend _NOT_ blindly calling normalized(u), since this will crash if |u| is 0.
            //       Instead, only update R_curr if (!IS_ZERO(squaredNorm(u)))
        } else if (input.mouse_left_released) {
            // TODO (see reading)
        }

        mat4 R_arcball = R_curr * R_prev; // you don't change have to this line


        // begin don't change this code
        { // NOTE I am disabling built-in camera rotation and 2D offsets
            camera_move(&camera);
            camera.theta = 0;
            camera.phi = 0;
            camera._o_x = 0;
            camera._o_y = 0;
        }
        mat4 P = camera_get_P(&camera);
        mat4 V = camera_get_V(&camera);
        if (mesh == 0) {
            fancy_draw(P, V, R_arcball, fancy_bunny, monokai.blue);
        } else {
            basic_draw(TRIANGLE_MESH, P * V * R_arcball, basic_box, monokai.green, 3, false, LERP(.8, monokai.green, monokai.white));
        }
        { // bespoke widget to visualize the arcball
            gl_begin(LINE_LOOP); {
                gl_color(monokai.orange);
                gl_PV(Identity4x4); // ! we are supplying vertices _in NDC_ (how cool is that?)
                for (double theta = 0; theta < 2 * PI - TINY; theta += (2 * PI) / 64) {
                    gl_vertex(e_theta(theta));
                }
            } gl_end();
        }
        // end   don't change this code

    }

}




int main() {
    init();
    hw4a();
    hw4b();
    hw4c();
    return 0;
}

