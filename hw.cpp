// #define COW_CRASH_ON_FLOATING_POINT_EXCEPTIONS
#define _CRT_SECURE_NO_WARNINGS
#include "snail.cpp"
#include "cow.cpp"
#include "_cow_supplement.cpp"

void hw6a() {
    init();

    // texture
    int side_length_in_pixels = 16;
    int nrChannels = 4; // r0 g0 b0 a0 r1 g1 b1 a1 ...
    unsigned char *data = (unsigned char *) malloc(side_length_in_pixels * side_length_in_pixels * nrChannels);
    memset(data, 255, side_length_in_pixels * side_length_in_pixels * nrChannels);

    fancy_texture_create("my_custom_texture", side_length_in_pixels, side_length_in_pixels, nrChannels, data);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST); // HACK so we can *see* the pixels
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST); // ""                             

    // mesh
    int3 triangle_indices[] = {
        { 0, 1, 2 },
        { 0, 2, 3}
    };
    vec3 vertex_positions[] = {
        { -1, -1, 0 },
        {  1, -1, 0 },
        {  1,  1, 0 },
        { -1,  1, 0 },
    };
    vec2 vertex_texCoords[] = {
        { 0, 0 },
        { 1, 0 },
        { 1, 1 },
        { 0, 1 }
    };

    Camera3D orbit = { 3, RAD(60) };
    double t = 0; FORNOW_UNUSED(t);
    while (begin_frame()) {
        t += 1. / 60;
        camera_move(&orbit);
        mat4 P = camera_get_P(&orbit);
        mat4 V = camera_get_V(&orbit);
        { // update texture
            // TODO match this reference  https://www.shadertoy.com/new
            // NOTE your solution must match the reference very closely
            fancy_texture_update("my_custom_texture", side_length_in_pixels, side_length_in_pixels, nrChannels, data);
        }
        fancy_draw(P, V, Identity4x4,
                2, triangle_indices,
                4, vertex_positions,
                NULL, NULL, {},
                vertex_texCoords, "my_custom_texture");
    }
}

void hw6b() {
    init();
    FPSCamera human = { V3(0, 10, 20), RAD(60), 0, 0 };
    while (begin_frame()) {
        fps_camera_move(&human);
        mat4 C = fps_camera_get_C(&human);
        mat4 P = tform_get_P_perspective(human.angle_of_view);
        mat4 V = inverse(C);
        mat4 PV = P * V; FORNOW_UNUSED(PV);
        { // pointer lock
            if (widget_active_widget_ID == 0 && input.mouse_left_pressed) {
                glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
            } else if (input.key_pressed[KEY_ESCAPE]) {
                glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
            }
        }
        { // // TODO skybox
            glDisable(GL_DEPTH_TEST); {
                int3 triangle_indices[] = { { 0, 1, 2 }, { 0, 2, 3} };
                vec2 vertex_texCoords[] = { { 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 } };

                // TODO set M to put the skybox on our head (it should _not_ rotate with us)
                mat4 M = Identity4x4;

                // NOTE the back of the skybox should be (initially) in _front_ of us
                //      i.e., i have set up back_vertex_positions correctly          
                // NOTE if you look down you can already see it                      
                vec3 back_vertex_positions[] = {
                    { -1, -1, -1 },
                    {  1, -1, -1 },
                    {  1,  1, -1 },
                    { -1,  1, -1 },
                };
                fancy_draw(P, V, M,
                        2, triangle_indices, 4, back_vertex_positions,
                        NULL, NULL, {},
                        vertex_texCoords, "data_texture_skybox_back.jpg");

                // TODO draw the front, left, right, bottom, and top

            } glEnable(GL_DEPTH_TEST);
        }
        { // TODO draw some small islands nearby
            // HINT if want something quick and dirty, meshlib.basic_tet will work
        }
        { // TODO draw the ocean
            // HINT i faked it with a single transparent quad, but follow your heart
        }
    }
}

void hw6c_draw_textured_square(mat4 P, mat4 V, mat4 M, char *texture_filename) {
    // please ignore; this function is hack-a-saurus rex
    static FancyTriangleMesh3D square;
    if (!square.num_vertices) {
        square = meshlib.fancy_square;
        square.vertex_normals = NULL;
    }
    fancy_draw(P, V, M,
            square.num_triangles, square.triangle_indices, square.num_vertices, square.vertex_positions,
            NULL, NULL, {},
            square.vertex_texCoords, texture_filename);
};

void hw6c() {
    init();

    struct {
        bool draw_depth_buffer = false;
        bool draw_scene_3D = true;
        bool draw_frustum = true;
        bool transparent_film_plane = false;
        bool draw_film_plane = true;
        bool draw_cube_at_observer = false;
        bool bunny_stress_test = false;
    } tweaks;

    struct {
        #define Q(deg, is_point) (2 * V3(cos(RAD(deg)), sin(RAD(deg)), (is_point) ? .5 : -.5))
        vec3 cycle_vertex_positions[9] = { Q(0, false), Q(20, false), Q(200, true), Q(120, false), Q(140, false), Q(320, true), Q(240, false), Q(260, false), Q(100, true),  };
        #undef Q
        vec3 cycle_vertex_colors[9] = { monokai.brown, monokai.brown, monokai.brown, monokai.yellow, monokai.yellow, monokai.yellow, monokai.purple, monokai.purple, monokai.purple, };
        BasicTriangleMesh3D cycle = { NELEMS(cycle_vertex_positions), cycle_vertex_positions, cycle_vertex_colors };

        vec3 tilt_vertex_positions[6] = { V3(-1, -.2, 1), V3(1, -.2, 1), V3(1, .2, -1), V3(-1, -.2, 1), V3(1, .2, -1), V3(-1, .2, -1), };
        vec3 tilt_vertex_colors[6] = { V3(1, 1, 1), V3(1, 1, 1), V3(0, 0, 1), V3(1, 1, 1), V3(0, 0, 1), V3(0, 0, 1) };
        BasicTriangleMesh3D tilt = { NELEMS(tilt_vertex_positions), tilt_vertex_positions, tilt_vertex_colors };

        vec3 clip2_vertex_positions[6] = { V3(-1, -1, 5), V3(1, -1, 5), V3(0, -1, 0) };
        vec3 clip2_vertex_colors[6] = { V3(1, 0, 0), V3(0, 1, 0), V3(0, 0, 1) };
        BasicTriangleMesh3D clip2 = { NELEMS(clip2_vertex_positions), clip2_vertex_positions, clip2_vertex_colors };

        vec3 clip1_vertex_positions[6] = { V3(-1, -1, 0), V3(1, -1, 0), V3(0, -1, 5) };
        vec3 clip1_vertex_colors[6] = { V3(1, 0, 0), V3(0, 1, 0), V3(0, 0, 1) };
        BasicTriangleMesh3D clip1 = { NELEMS(clip1_vertex_positions), clip1_vertex_positions, clip1_vertex_colors };
    } meshes;
    BasicTriangleMesh3D *examples[] = { &meshes.cycle, &meshes.tilt, &meshes.clip2, &meshes.clip1 };
    char *example_names[] = { "cycle", "tilt", "clip2", "clip1" };
    BasicTriangleMesh3D bunny = load_basic_mesh("data_basic_bunny", true);

    int side_length_in_pixels = 64;

    Texture depth_buffer = { "depth_buffer", side_length_in_pixels, side_length_in_pixels, 1, (unsigned char *) malloc(side_length_in_pixels * side_length_in_pixels * 1) };
    Texture color_buffer = { "color_buffer", side_length_in_pixels, side_length_in_pixels, 4, (unsigned char *) malloc(side_length_in_pixels * side_length_in_pixels * 4) };

    fancy_texture_create(&depth_buffer); glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST); glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    fancy_texture_create(&color_buffer); glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST); glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

    Camera3D observer = { 7, RAD(60), RAD(30), RAD(-15), -2 };
    Camera3D renderer = { 3, RAD(60) };
    double renderer_n = -1; // near clip plane
    double renderer_f = -5; // far clip plane 
    double renderer_distance_to_film_plane = -renderer_n; // for visualization only (doesn't impact rendering)

    int mesh_index = 0;
    while (begin_frame()) {
        BasicTriangleMesh3D *mesh = examples[mesh_index]; {
            imgui_slider("mesh_index", &mesh_index, 0, NELEMS(examples) - 1, 'j', 'k', true);
            if (!tweaks.bunny_stress_test) {
                _imgui_printf("example: %s", example_names[mesh_index]);
            } else {
                mesh = &bunny;
                _imgui_printf("i'm starting to worry our implementation wasn't very efficient");
            }
        }

        { // tweaks
            imgui_checkbox("draw_depth_buffer", &tweaks.draw_depth_buffer, 'z');
            imgui_checkbox("draw_film_plane", &tweaks.draw_film_plane, 'x');
            imgui_checkbox("draw_frustum", &tweaks.draw_frustum, 'f');
            imgui_checkbox("draw_scene_3D", &tweaks.draw_scene_3D, 's');
            imgui_checkbox("draw_cube_at_observer", &tweaks.draw_cube_at_observer, 'c');
            imgui_checkbox("transparent_film_plane", &tweaks.transparent_film_plane, 't');
            imgui_checkbox("bunny_stress_test", &tweaks.bunny_stress_test);

            imgui_slider("renderer_distance_to_film_plane", &renderer_distance_to_film_plane, -renderer_n, -renderer_f);
            imgui_slider("renderer.distance_to_origin", &renderer.distance_to_origin, 1, 10);
            imgui_slider("renderer.theta", &renderer.theta, RAD(-180), RAD(180), true);
            imgui_slider("renderer.phi", &renderer.phi, RAD(-90), RAD(90), true);
            imgui_slider("renderer.angle_of_view", &renderer.angle_of_view, RAD(1), RAD(178), true);
            imgui_slider("renderer_n", &renderer_n, -2, -.1);
            imgui_slider("renderer_f", &renderer_f, -20, renderer_n);
            renderer_distance_to_film_plane = CLAMP(renderer_distance_to_film_plane, -renderer_n, -renderer_f);
        }

        mat4 C_rasterizer = camera_get_C(&renderer);

        { // // // TODO rasterize
            mat4 P_rasterizer = tform_get_P_perspective(renderer.angle_of_view, renderer_n, renderer_f, 1);
            mat4 V_rasterizer = inverse(C_rasterizer);
            { // init buffers
                memset(depth_buffer.data, 255, side_length_in_pixels * side_length_in_pixels);
                if (tweaks.transparent_film_plane) {
                    memset(color_buffer.data, 0, side_length_in_pixels * side_length_in_pixels * 4);
                } else {
                    memset(color_buffer.data, 255, side_length_in_pixels * side_length_in_pixels * 4);
                    for (int i = 0; i < side_length_in_pixels * side_length_in_pixels; ++i) {
                        color_buffer.data[4 * i + 3] = 190;
                    }
                }
            }
            for (int triangle_i = 0; triangle_i < mesh->num_vertices / 3; ++triangle_i) {
                vec3 a_world = mesh->vertex_positions[3 * triangle_i + 0];
                vec3 b_world = mesh->vertex_positions[3 * triangle_i + 1];
                vec3 c_world = mesh->vertex_positions[3 * triangle_i + 2];
                vec3 a_camera = transformPoint(V_rasterizer, a_world);
                vec3 b_camera = transformPoint(V_rasterizer, b_world);
                vec3 c_camera = transformPoint(V_rasterizer, c_world);
                vec2 a_NDC = transformPoint(P_rasterizer, a_camera).xy;
                vec2 b_NDC = transformPoint(P_rasterizer, b_camera).xy;
                vec2 c_NDC = transformPoint(P_rasterizer, c_camera).xy;
                vec3 fallback_color = V3(1, 0, 1); {
                    // TODOLATER fallback__color <- V3(.5, .5, .5) + .5 * face_normal (in world coords)
                    // NOTE the bunny should now look rather wonderful                                 
                }
                vec3 color_a = (mesh->vertex_colors != NULL) ? mesh->vertex_colors[3 * triangle_i + 0] : fallback_color;
                vec3 color_b = (mesh->vertex_colors != NULL) ? mesh->vertex_colors[3 * triangle_i + 1] : fallback_color;
                vec3 color_c = (mesh->vertex_colors != NULL) ? mesh->vertex_colors[3 * triangle_i + 2] : fallback_color;
                FORNOW_UNUSED(a_NDC); FORNOW_UNUSED(b_NDC); FORNOW_UNUSED(c_NDC);
                FORNOW_UNUSED(color_a); FORNOW_UNUSED(color_b); FORNOW_UNUSED(color_c);
                // // TODO rasterize                                                                   
                // // NOTE start this problem early! (it may be a doozy)                               
                // // HINT be careful with unsigned char's (don't, e.g., use a regular char by mistake)
                // // HINT my solution (to everything except clipping) is 40 additional lines          
                //                                                                                     
                // TODO iterate over all the pixels and color blue if inside the triangle              
                // NOTE the triangles should show up for the cycle and tilt examples                   
                // NOTE ignore the clip examples                                                       
                // HINT color_buffer.data is an unsigned char array storing r0 g0 b0 a0 r1 g1 b1 a1 ...
                // HINT a of 0 is fully transparent; a of 255 is fully opaque                          
                // HINT (r, g, b, a) <- (0, 0, 255, 255) is full blue                                  
                //                                                                                     
                // TODO color the pixel a more appropriate color using barycentric interpolation in NDC
                // NOTE the cycle example colors will be correct, but triangles will overlap wrong     
                // NOTE the colors of the tilt example will look (slightly) wrong                      
                //                                                                                     
                // TODO write and use depth_buffer.data to make the overlapping correct                
                // NOTE "Depth Buffer" and "Z Buffer" are different names for the same thing           
                // NOTE the cycle example is now complete :D                                           
                // HINT z buffer of 0 is at the near clip plane, z buffer of 255 is at far clip plane  
                // HINT to get z in camera coordinates, please see: lowk_persp_interp_techrep.pdf      
                //      HINT where he uses linear interpolation (between two vertices),                
                //           you can substitute barycentric interpolation (between three vertices)     
                //                                                                                     
                // TODO perspective-correct color interpolation                                        
                // NOTE the tilt example is now complete :D                                            
                // HINT again see this writeup: lowk_persp_interp_techrep.pdf                          
                //                                                                                     
                // TODO (if you want an extra challenge) implement clipping                            
                // NOTE the clip examples should hopefully no longer be broken?                        

            }
            fancy_texture_update(&depth_buffer);
            fancy_texture_update(&color_buffer);
        }

        { // observe
            camera_move(&observer);
            char *name = "hold to toggle C_observer <- C_rasterizer";

            static bool clicked;
            bool clicked_this_frame;
            bool released_this_frame = false;
            {
                clicked_this_frame = imgui_button(name, KEY_TAB);
                if (clicked_this_frame) {
                    clicked = true;
                }
                bool selected = (imgui.selected_widget_ID == (void *) name);
                if (clicked && !selected) {
                    clicked = false;
                    released_this_frame = true;
                }
            }
            {
                static Camera3D safe;
                static bool draw_cube_push;
                if (clicked_this_frame) {
                    memcpy(&safe, &observer, sizeof(Camera3D));
                    memcpy(&observer, &renderer, sizeof(Camera3D));
                    draw_cube_push = tweaks.draw_cube_at_observer;
                    tweaks.draw_cube_at_observer = false;
                }
                if (released_this_frame) {
                    memcpy(&observer, &safe, sizeof(Camera3D));
                    tweaks.draw_cube_at_observer = draw_cube_push;
                }
            }
            mat4 P_observer = camera_get_P(&observer);
            mat4 V_observer = camera_get_V(&observer);
            mat4 PV_observer = P_observer * V_observer;

            if (tweaks.draw_scene_3D) {
                basic_draw(TRIANGLES, PV_observer, *mesh, V3(1, 0, 1));
            }

            { // bespoke widget
                basic_draw(PV_observer * C_rasterizer, meshlib.basic_axes);

                double n = renderer_n;
                double f = renderer_f;
                double hangle = renderer.angle_of_view / 2;
                double r_n = -n * tan(hangle);
                double r_f = -f * tan(hangle);
                if (tweaks.draw_frustum) {
                    vec3 frustum_vertex_positions[] = {{r_n,r_n,n},{-r_n,r_n,n},{-r_n,-r_n,n},{r_n,-r_n,n},{r_n,r_n,n},{r_f,r_f,f},{-r_f,r_f,f},{-r_f,-r_f,f},{r_f,-r_f,f},{r_f,r_f,f},{r_f,-r_f,f},{r_n,-r_n,n},{r_f,-r_f,f},{-r_f,-r_f,f},{-r_n,-r_n,n},{-r_f,-r_f,f},{-r_f,r_f,f},{-r_n,r_n,n},{-r_f,r_f,f},{r_f,r_f,f},{r_n,r_n,n},};
                    basic_draw(LINE_STRIP, PV_observer * C_rasterizer, NELEMS(frustum_vertex_positions), frustum_vertex_positions, monokai.white, 3);
                }
                if (tweaks.draw_film_plane) {
                    double s = LERP(INVERSE_LERP(renderer_distance_to_film_plane, -n, -f), r_n, r_f);
                    mat4 M = C_rasterizer * Translation(0, 0, -renderer_distance_to_film_plane) * Scaling(s);
                    hw6c_draw_textured_square(P_observer, V_observer, M, (!tweaks.draw_depth_buffer) ? color_buffer.filename : depth_buffer.filename);
                    { // outline
                        vec3 tmp[] = { { -1, -1, 0 }, { -1,  1, 0 }, {  1,  1, 0 }, {  1, -1, 0 }, };
                        basic_draw(LINE_LOOP, P_observer * V_observer * M, 4, tmp);
                    }
                }
                if (tweaks.draw_cube_at_observer) {
                    basic_draw(PV_observer * C_rasterizer * Scaling(.25), meshlib.basic_box, .5 * monokai.gray);
                }
            }
        }
    }
}

int main() {
    hw6a();
    hw6b();
    hw6c();
    return 0;
}

