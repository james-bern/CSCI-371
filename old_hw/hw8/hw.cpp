#define _CRT_SECURE_NO_WARNINGS
#include "snail.cpp"
#include "cow.cpp"
#include "_cow_supplement.cpp"

// BEGIN fine to ignore                                                         
void hw8a_draw_textured_square(mat4 P, mat4 V, mat4 M, char *texture_filename) {
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

struct {
    #define Q(deg) (2 * V3(cos(RAD(deg)), sin(RAD(deg)), -2))
    vec3 trivial_vertex_positions[3] = { Q(0), Q(120), Q(240)};
    #undef Q
    vec3 trivial_vertex_colors[3] = { V3(1, 0, 0), V3(0, 1, 0), V3(0, 0, 1) };
    BasicTriangleMesh3D trivial = { NELEMS(trivial_vertex_positions), trivial_vertex_positions, trivial_vertex_colors };

    #define Q(deg, is_point) (2 * V3(cos(RAD(deg)), sin(RAD(deg)), (is_point) ? .5 : -.5))
    vec3 cycle_vertex_positions[12] = { Q(0, false), Q(20, false), Q(200, true), Q(120, false), Q(140, false), Q(320, true), Q(240, false), Q(260, false), Q(100, true), V3(5 * cos(RAD(240)), -1.5, 5 * sin(RAD(240))), V3(5 * cos(RAD(120)), -1.5, 5 * sin(RAD(120))), V3(5 * cos(RAD(0)), -1.5, 5 * sin(RAD(0))) };
    #undef Q
    vec3 cycle_vertex_colors[12] = { monokai.brown, monokai.brown, monokai.brown, monokai.yellow, monokai.yellow, monokai.yellow, monokai.purple, monokai.purple, monokai.purple, monokai.red, monokai.red, monokai.red };
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
} hw8a_meshes;

struct {
    bool draw_rays;
    bool draw_scene_3D = true;
    bool fully_transparent_film_plane_bg;
    bool draw_film_plane = true;
    bool draw_cube_at_observer;
    bool bunny_stress_test;
    double renderer_distance_to_film_plane = 1; // for visualization only (doesn't impact rendering)
} hw8a_tweaks;
// END fine to ignore                                                           




// mesh    -- the current scene
// light_p -- the position of the point light in world coordinates
BasicTriangleMesh3D *mesh;               
vec3 light_p = V3(0, 2.5, 0);

// S                  -- the film plane side length in pixels
// color_buffer.data  -- an unsigned char array in r_0 g_0 b_0 a_0 r_1 ... format
//                       where, e.g., g_k is the green component of the k-th pixel
// hw8a_set_pixel(...) -- write color to pixel (i, j); !! does clamping for you
int S = 8;
Texture color_buffer;
void hw8a_set_pixel(int i, int j, vec3 color) {
    color_buffer.data[4 * (j * S + i) + 0] = (unsigned char)(255 * CLAMP(color.r, 0, 1));
    color_buffer.data[4 * (j * S + i) + 1] = (unsigned char)(255 * CLAMP(color.g, 0, 1));
    color_buffer.data[4 * (j * S + i) + 2] = (unsigned char)(255 * CLAMP(color.b, 0, 1));
    color_buffer.data[4 * (j * S + i) + 3] = 255;
}

void hw8a() {
    init();

    // renderer -- the camera doing the actual rendering (raycasting)
    // observer -- the camera that's observing _everything_ (scene and ray-casting camera)
    Camera3D renderer = { 4, RAD(60) };
    Camera3D observer = { 6.5, RAD(60), RAD(30), RAD(-15), -2 };
    { // (fine to ignore)
        color_buffer= { "color_buffer", S, S, 4, (unsigned char *) malloc(S * S * 4) };
        fancy_texture_create(&color_buffer); glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST); glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    }

    while (begin_frame()) {
        { // tweaks (fine to ignore)
            imgui_checkbox("draw_rays", &hw8a_tweaks.draw_rays, 'z');
            imgui_checkbox("draw_film_plane", &hw8a_tweaks.draw_film_plane, 'x');
            imgui_checkbox("draw_scene_3D", &hw8a_tweaks.draw_scene_3D, 's');
            imgui_checkbox("draw_cube_at_observer", &hw8a_tweaks.draw_cube_at_observer, 'c');
            imgui_checkbox("fully_transparent_film_plane_bg", &hw8a_tweaks.fully_transparent_film_plane_bg, 't');
            imgui_checkbox("bunny_stress_test", &hw8a_tweaks.bunny_stress_test);

            imgui_slider("hw8a_tweaks.renderer_distance_to_film_plane", &hw8a_tweaks.renderer_distance_to_film_plane, 0, 5);
            imgui_slider("renderer.distance_to_origin", &renderer.distance_to_origin, 1, 10);
            imgui_slider("renderer.theta", &renderer.theta, RAD(-30), RAD(30), true);
            imgui_slider("renderer.phi", &renderer.phi, RAD(-30), RAD(30), true);
            imgui_slider("renderer.angle_of_view", &renderer.angle_of_view, RAD(1), RAD(178), true);
        }

        // *_renderer                   -- matrices for the renderer
        // *_observere                  -- matrices for the observer
        // film_plane_side_length_world -- the length of the film plane in world coordinates
        mat4 C_renderer, P_renderer, V_renderer; 
        mat4 P_observer, V_observer, PV_observer;
        { // (fine to ignore)
            { // mesh
                static BasicTriangleMesh3D *examples[] = { &hw8a_meshes.trivial, &hw8a_meshes.cycle, &hw8a_meshes.tilt, &hw8a_meshes.clip2, &hw8a_meshes.clip1 };
                static char *example_names[] = { "trivial", "cycle", "tilt", "clip2", "clip1" };
                static BasicTriangleMesh3D bunny = load_basic_mesh("data_basic_bunny", true);
                static int mesh_index = 0;
                imgui_slider("mesh_index", &mesh_index, 0, NELEMS(examples) - 1, 'j', 'k', true);
                mesh = examples[mesh_index];
                if (!hw8a_tweaks.bunny_stress_test) {
                    _imgui_printf("example: %s", example_names[mesh_index]);
                } else {
                    mesh = &bunny;
                    _imgui_printf("i'm starting to worry our implementation wasn't very efficient");
                }
            }

            { // C_*, P_*, PV_*
                { // reset
                    static Camera3D _renderer_0 = renderer;
                    if (imgui_button("reset C_renderer", 'r')) {
                        renderer = _renderer_0;
                    }
                }
                camera_move(&observer);
                C_renderer = camera_get_C(&renderer);
                P_renderer = tform_get_P_perspective(renderer.angle_of_view); // aspect <- 1
                V_renderer = inverse(C_renderer);
                { // C_observer <- C_renderer
                    static bool clicked;
                    bool clicked_this_frame;
                    bool selected;
                    bool released_this_frame = false;
                    {
                        char *name = "hold to toggle C_observer <- C_renderer";
                        clicked_this_frame = imgui_button(name, KEY_TAB);
                        if (clicked_this_frame) {
                            clicked = true;
                        }
                        selected = (imgui.selected_widget_ID == (void *) name);
                        if (clicked && !selected) {
                            clicked = false;
                            released_this_frame = true;
                        }
                    }

                    { // memcpy's
                        static Camera3D safe;
                        static bool draw_cube_push;
                        if (clicked_this_frame) {
                            memcpy(&safe, &observer, sizeof(Camera3D));
                            draw_cube_push = hw8a_tweaks.draw_cube_at_observer;
                            hw8a_tweaks.draw_cube_at_observer = false;
                        }
                        if (selected) {
                            memcpy(&observer, &renderer, sizeof(Camera3D));
                        }
                        if (released_this_frame) {
                            memcpy(&observer, &safe, sizeof(Camera3D));
                            hw8a_tweaks.draw_cube_at_observer = draw_cube_push;
                        }
                    }
                }
                P_observer = camera_get_P(&observer);
                V_observer = camera_get_V(&observer);
                PV_observer = P_observer * V_observer;
            }
        }


        { // render (your work here!)
            { // prep (fine to ignore)
                { // clear color_buffer
                    if (hw8a_tweaks.fully_transparent_film_plane_bg) {
                        memset(color_buffer.data, 0, S * S * 4);
                    } else {
                        memset(color_buffer.data, 255, S * S * 4);
                        int pixelIndex = 0;
                        for (int j = 0; j < S; ++j) {
                            for (int i = 0; i < S; ++i) {
                                color_buffer.data[4 * pixelIndex++ + 3] = (((j + i) % 2) == 0) ? 200 : 220;
                            }
                        }
                    }
                }

                { // light_p
                    jank_widget_translate3D(PV_observer, 1, &light_p);
                    basic_draw(POINTS, PV_observer, 1, &light_p, monokai.white);
                }
            }

            { // render (your work here!)
                // renderer.angle_of_view -- _full_ (not half) angle of view of the renderer
                // *_renderer             -- the axes and origin of the renderer
                // NOTE: o_renderer is where rays originate from
                // NOTE: -z_renderer points from o_renderer to the center of the film plane
                vec3 x_renderer, y_renderer, z_renderer, o_renderer;
                {
                    camera_get_coordinate_system(&renderer, NULL, x_renderer.data, y_renderer.data, z_renderer.data, o_renderer.data);
                }

                { // write to color_buffer.data (your work here!)

                    // gl_* will be useful for debugging direction (dir)
                    gl_begin(LINES); gl_PV(PV_observer); {


                        for (int i = 0; i < S; ++i) {
                            for (int j = 0; j < S; ++j) {

                                double theta = renderer.angle_of_view;
                                vec3 o = o_renderer;
                                vec3 dir = V3(i - double(S) / 2, j - double(S) / 2, -(double(S) / 2) / tan(theta / 2)); // student answer from board
                                if (hw8a_tweaks.draw_rays) {
                                    gl_color(monokai.red);
                                    gl_vertex(o);
                                    gl_vertex(o + dir);
                                }

                                // TODO intersect with the scene
                                // for (...) {
                                //     if (...) {
                                //         hw8a_set_pixel(i, j, color);
                                //     }
                                // }
                            }
                        }

                    } gl_end();

                }

                { // send updated texture to the GPU (fine to ignore)
                    fancy_texture_update(&color_buffer);
                }
            }
        }

        { // observe (fine to ignore)
            if (hw8a_tweaks.draw_scene_3D) {
                basic_draw(TRIANGLES, PV_observer, *mesh, V3(1, 0, 1));
            }
            { // bespoke widget
                basic_draw(PV_observer * C_renderer, meshlib.basic_axes);

                if (hw8a_tweaks.draw_film_plane) {
                    double film_plane_side_length_world = 2 * hw8a_tweaks.renderer_distance_to_film_plane * tan(renderer.angle_of_view / 2);
                    mat4 M = C_renderer * Translation(0, 0, -hw8a_tweaks.renderer_distance_to_film_plane) * Scaling(film_plane_side_length_world / 2);
                    hw8a_draw_textured_square(P_observer, V_observer, M, color_buffer.filename);
                    { // outline
                        vec3 tmp[] = { { -1, -1, 0 }, { -1,  1, 0 }, {  1,  1, 0 }, {  1, -1, 0 }, };
                        basic_draw(LINE_LOOP, P_observer * V_observer * M, 4, tmp);
                    }
                }
                if (hw8a_tweaks.draw_cube_at_observer) {
                    basic_draw(PV_observer * C_renderer * Scaling(.2), meshlib.basic_box, .5 * monokai.gray);
                }
            }
        }
    }
}

char *hw8b_frag = R""""(
    #version 330 core

    vec3 plasma(float t) {
        const vec3 c0 = vec3(0.05873234392399702, 0.02333670892565664, 0.5433401826748754);
        const vec3 c1 = vec3(2.176514634195958, 0.2383834171260182, 0.7539604599784036);
        const vec3 c2 = vec3(-2.689460476458034, -7.455851135738909, 3.110799939717086);
        const vec3 c3 = vec3(6.130348345893603, 42.3461881477227, -28.51885465332158);
        const vec3 c4 = vec3(-11.10743619062271, -82.66631109428045, 60.13984767418263);
        const vec3 c5 = vec3(10.02306557647065, 71.41361770095349, -54.07218655560067);
        const vec3 c6 = vec3(-3.658713842777788, -22.93153465461149, 18.19190778539828);
        return c0+t*(c1+t*(c2+t*(c3+t*(c4+t*(c5+t*c6)))));
    }

    // https://iquilezles.org/articles/distfunctions/
    float dot2(vec2 v) { return dot(v,v); }
    float dot2(vec3 v) { return dot(v,v); }
    float ndot(vec2 a, vec2 b) { return a.x*b.x - a.y*b.y; }

    // for computing ray directions
    uniform vec3 x_renderer;
    uniform vec3 y_renderer;
    uniform vec3 z_renderer;
    uniform vec3 o_renderer;
    uniform float renderer_angle_of_view;
    uniform vec2 iResolution;

    uniform float time; // for time-varying distance fields

    out vec4 fragColor;

    // begin https://iquilezles.org/articles/distfunctions/
    float sdSphere(vec3 p, float r) {
        return length(p) - r;
    }
    float sdTorus(vec3 p, vec2 t) {
        vec2 q = vec2(length(p.xz)-t.x,p.y);
        return length(q)-t.y;
    }
    // end

    vec3 march(vec3 o, vec3 dir) {
        // https://michaelwalczyk.com/blog-ray-marching.html

        const int MAX_STEPS = 64;
        const float HIT_TOLERANCE = 0.001;
        const float MAX_MARCH_DISTANCE = 100.0;

        // p   -- current position along ray
        // o   -- camera origin             
        // t   -- distance marched along ray
        // dir -- camera direction          
        // f   -- distance to surface       

        float t = 0.0;
        int step = 0;
        while (step++ < MAX_STEPS && t < MAX_MARCH_DISTANCE) {
            // get current position of ray's head
            vec3 p = o + t * dir;

            // compute distance to implicit surface
            float f = MAX_MARCH_DISTANCE; {
                {
                    vec3 sphere_position = vec3(0.0, 0.0, 0.0);
                    float sphere_radius = 1.0;
                    float distance_to_sphere = sdSphere(p - sphere_position, sphere_radius);
                    f = min(f, distance_to_sphere);
                }
                {
                    // f = min(f, sdTorus(p - vec3(0.0, 1.0 * sin(time), 0.0), vec2(1.0, 0.25)));
                    vec3 torus_position = vec3(0.0, sin(time), 0.0);
                    float torus_major_radius = 1.0;
                    float torus_minor_radius = 0.25;
                    vec2 torus_radii = vec2(torus_major_radius, torus_minor_radius);
                    float distance_to_torus = sdTorus(p - torus_position, torus_radii);
                    f = min(f, distance_to_torus);
                }
            }

            if (f < HIT_TOLERANCE) { // hit!
                return plasma(.5 + .5 * cos(t));
            }

            // NOTE if you're getting weird "overstepping" artifacts
            // (weird missing slices in the geometry)               
            // a (hacky) solution is to replace t += f; with e.g.   
            // t += min(f, .5);
            t += f;
        }
        return vec3(0.0);
    }

    void main() {
        vec3 o = o_renderer;
        vec3 dir; {
            // NOTE assume unit distance to film plane
            vec2 ds; { // [-R, R]
                float theta = renderer_angle_of_view / 2;
                float _R = tan(theta);
                ds = gl_FragCoord.xy;
                ds -= vec2(iResolution.x / 2, iResolution.y / 2);
                ds *= _R * 2. / iResolution.y;
            }
            // vec3 p_world = o_renderer - z_renderer + dx * x_renderer + dy * y_renderer;
            // dir = p_world - o_renderer;
            dir = -z_renderer + ds.x * x_renderer + ds.y * y_renderer;
        }
        vec3 col = march(o, dir);
        fragColor = vec4(col, 1);
    }
)"""";

void hw8b() {
    init();

    // mesh
    int num_vertices = 4;
    int num_triangles = 2;
    vec3 vertex_positions[] = { { -1, -1, 0 }, { 1, -1, 0 }, { 1, 1, 0 }, { -1, 1, 0 } };
    int3 triangle_indices[] = { { 0, 1, 2 }, { 0, 2, 3 } };

    // shaders
    char *vert = R""""(
        #version 330 core
        layout (location = 0) in vec3 _p_model;
        void main() {
            gl_Position = vec4(_p_model, 1);
        }
    )"""";

    int shader_program = shader_build_program(vert, hw8b_frag);
    ASSERT(shader_program);

    // misc opengl
    GLuint VAO, VBO, EBO; {
        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);
        glGenBuffers(1, &EBO);
    }

    // misc cow
    bool playing = false;
    double time = 0;
    Camera3D renderer = { 8, RAD(45) };

    while (begin_frame()) {
        camera_move(&renderer);

        { // imgui
            { // reset
                static Camera3D _camera_0 = renderer;
                if (imgui_button("reset", 'r')) {
                    renderer = _camera_0;
                    time = 0;
                }
            }
            imgui_checkbox("playing", &playing, 'p');
        }

        { // draw
            glBindVertexArray(VAO);
            glBindBuffer(GL_ARRAY_BUFFER, VBO);
            glBufferData(GL_ARRAY_BUFFER, num_vertices * 3 * sizeof(double), vertex_positions, GL_DYNAMIC_DRAW);
            glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 0, NULL);
            glEnableVertexAttribArray(0);

            glUseProgram(shader_program);

            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3 * num_triangles * sizeof(int), triangle_indices, GL_DYNAMIC_DRAW);

            shader_set_uniform_vec2(shader_program, "iResolution", window_get_dimensions_in_pixels());
            {
                shader_set_uniform_double(shader_program, "time", time);

                shader_set_uniform_double(shader_program, "renderer_angle_of_view", renderer.angle_of_view);
                vec3 x_renderer, y_renderer, z_renderer, o_renderer; {
                    camera_get_coordinate_system(&renderer, NULL, x_renderer.data, y_renderer.data, z_renderer.data, o_renderer.data);
                }
                shader_set_uniform_vec3(shader_program, "x_renderer", x_renderer);
                shader_set_uniform_vec3(shader_program, "y_renderer", y_renderer);
                shader_set_uniform_vec3(shader_program, "z_renderer", z_renderer);
                shader_set_uniform_vec3(shader_program, "o_renderer", o_renderer);
            }


            glDrawElements(GL_TRIANGLES, 3 * num_triangles, GL_UNSIGNED_INT, NULL);
        }

        if (playing || input.key_pressed['.']) {
            time += .0167;
        }
    }
}



int main() {
    hw8a();
    hw8b();
    return 0;
}

