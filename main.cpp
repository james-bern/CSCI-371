// #define COW_PATCH_FRAMERATE
// #define COW_PATCH_FRAMERATE_SLEEP
#include "include.cpp"

////////////////////////////////////////////////////////////////////////////////
// software rasterizer /////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void rasterize(
        mat4 P,
        mat4 V,
        mat4 M,
        IndexedTriangleMesh3D *mesh,
        Texture *color_buffer,
        Texture *depth_buffer,
        real z_near,
        real z_far) {

    // TODO clear the depth_buffer to 1.0 (max depth)
    // TODO clear the color_buffer to whatever you like; i like rgb = V3(1.0), a = 0.5
    int side_length_in_pixels = color_buffer->width;
    for (int i = 0; i < side_length_in_pixels; ++i) {
        for (int j = 0; j < side_length_in_pixels; ++j) {
            texture_set_pixel(depth_buffer, i, j, 1.0);
            texture_set_pixel(color_buffer, i, j, V3(1.0, 1.0, 1.0), .5);
        }
    }

    // TODO: iterate over all triangles and rasterize
    // TODO: if the mesh doesn't have vertex_colors specified
    //       (i.e. if (mesh->vertex_colors == NULL) { ... })
    //       please use (V3(.5) + .5 * vertex_normal) as the vertex color
    for (int triangle_i = 0; triangle_i < mesh->num_triangles; ++triangle_i) {
        vec3 a_world = transformPoint(M, mesh->vertex_positions[mesh->triangle_indices[triangle_i][0]]);
        vec3 b_world = transformPoint(M, mesh->vertex_positions[mesh->triangle_indices[triangle_i][1]]);
        vec3 c_world = transformPoint(M, mesh->vertex_positions[mesh->triangle_indices[triangle_i][2]]);
        vec3 a_camera = transformPoint(V, a_world);
        vec3 b_camera = transformPoint(V, b_world);
        vec3 c_camera = transformPoint(V, c_world);
        vec3 a_NDC = transformPoint(P, a_camera);
        vec3 b_NDC = transformPoint(P, b_camera);
        vec3 c_NDC = transformPoint(P, c_camera);
        vec3 color_a, color_b, color_c;
        if (mesh->vertex_colors != NULL) {
            color_a = mesh->vertex_colors[mesh->triangle_indices[triangle_i][0]];
            color_b = mesh->vertex_colors[mesh->triangle_indices[triangle_i][1]];
            color_c = mesh->vertex_colors[mesh->triangle_indices[triangle_i][2]];
        } else {
            color_a = V3(.5) + .5 * mesh->vertex_normals[mesh->triangle_indices[triangle_i][0]];
            color_b = V3(.5) + .5 * mesh->vertex_normals[mesh->triangle_indices[triangle_i][1]];
            color_c = V3(.5) + .5 * mesh->vertex_normals[mesh->triangle_indices[triangle_i][2]];
        }
        mat3 S = {
            a_NDC.x, b_NDC.x, c_NDC.x,
            a_NDC.y, b_NDC.y, c_NDC.y,
            1,   1,   1
        };
        for (int i = 0; i < side_length_in_pixels; ++i) {
            real p_y_NDC = -1 + 2 * real(i) / (side_length_in_pixels - 1);
            if (p_y_NDC > MAX(MAX(a_NDC.y, b_NDC.y), c_NDC.y)) { continue; }
            if (p_y_NDC < MIN(MIN(a_NDC.y, b_NDC.y), c_NDC.y)) { continue; }
            for (int j = 0; j < side_length_in_pixels; ++j) {
                real p_x_NDC = -1 + 2 * real(j) / (side_length_in_pixels - 1);
                vec3 p_NDC = { p_x_NDC, p_y_NDC, 1 };
                if (p_x_NDC > MAX(MAX(a_NDC.x, b_NDC.x), c_NDC.x)) { continue; }
                if (p_x_NDC < MIN(MIN(a_NDC.x, b_NDC.x), c_NDC.x)) { continue; }
                vec3 weights_NDC = inverse(S) * p_NDC;
                if (weights_NDC[0] > 0 && weights_NDC[1] > 0 && weights_NDC[2] > 0) {
                    real z_camera = weights_NDC[0] * a_camera.z + weights_NDC[1] * b_camera.z + weights_NDC[2] * c_camera.z;
                    // real z_camera = 1 / (weights_NDC[0] * (1 / a_camera.z) + weights_NDC[1] * (1 / b_camera.z) + weights_NDC[2] * (1 / c_camera.z));
                    if (z_far < z_camera && z_camera < z_near) {
                        real z = INVERSE_LERP(1 / z_camera, 1 / z_near, 1 / z_far);
                        real curr; texture_get_pixel(depth_buffer, i, j, &curr);
                        if (z < curr) {
                            texture_set_pixel(depth_buffer, i, j, z);
                            vec3 color = weights_NDC[0] * color_a + weights_NDC[1] * color_b + weights_NDC[2] * color_c;
                            // vec3 color = (weights_NDC[0] * (color_a / a_camera.z) + weights_NDC[1] * (color_b / b_camera.z) + weights_NDC[2] * (color_c / c_camera.z)) / (1 / z_camera);
                            texture_set_pixel(color_buffer, i, j, color, 1.0);
                        }
                    }
                }
            }
        }
    }
}

vec3 _example_vertex_positions[] = {
    { cos(RAD(000)), sin(RAD(000)), -1.0 },
    { cos(RAD(020)), sin(RAD(020)), -1.0 },
    { cos(RAD(200)), sin(RAD(200)),  1.0 },
    { cos(RAD(120)), sin(RAD(120)), -1.0 },
    { cos(RAD(140)), sin(RAD(140)), -1.0 },
    { cos(RAD(320)), sin(RAD(320)),  1.0 },
    { cos(RAD(240)), sin(RAD(240)), -1.0 },
    { cos(RAD(260)), sin(RAD(260)), -1.0 },
    { cos(RAD(100)), sin(RAD(100)),  1.0 },
};
vec3 _example_vertex_colors[] = {
    monokai.yellow,
    monokai.yellow,
    monokai.yellow,
    monokai.purple,
    monokai.purple,
    monokai.purple,
    monokai.brown,
    monokai.brown,
    monokai.brown,
};
int3 _example_triangle_indices[] = {
    { 0, 1, 2 },
    { 3, 4, 5 },
    { 6, 7, 8 },
};

void hw5a() {

    IndexedTriangleMesh3D mesh_example = {}; {
        mesh_example.num_vertices = 9;
        mesh_example.num_triangles = 3;
        mesh_example.vertex_positions = _example_vertex_positions;
        mesh_example.vertex_colors = _example_vertex_colors;
        mesh_example.triangle_indices = _example_triangle_indices;
    }

    int side_length_in_pixels = 64;
    Texture color_buffer = texture_create("color_buffer", side_length_in_pixels, side_length_in_pixels, 4);
    Texture depth_buffer = texture_create("depth_buffer", side_length_in_pixels, side_length_in_pixels, 1);

    Camera3D observer = { 7, RAD(30), RAD(30), RAD(-15), -2 };
    Camera3D renderer = { 3, RAD(30) };
    real renderer_n = -2; // near clip plane
    real renderer_f = -4; // far clip plane

    bool draw_bunny_instead = false;
    bool draw_depth_instead_of_color = false;
    bool _hide_camera_cube = false;

    real _renderer_distance_to_film_plane = -renderer_n; // for visualization only (doesn't impact rendering)
    bool playing = false;
    real time = 0.0;
    while (cow_begin_frame()) {
        if (draw_bunny_instead) {
            gui_checkbox("playing", &playing, 'p');
        }
        if (playing) { time += .0167; }
        gui_checkbox("draw_depth_instead_of_color", &draw_depth_instead_of_color, 'z');
        gui_checkbox("draw_bunny_instead", &draw_bunny_instead, 'b');
        gui_checkbox("_hide_camera_cube", &_hide_camera_cube, 'c');

        IndexedTriangleMesh3D *mesh = &mesh_example;
        mat4 M = globals.Identity;
        if (draw_bunny_instead) {
            mesh = &library.meshes.bunny;
            M = M4_RotationAboutYAxis(time) * M4_Scaling(.7);
        }

        { // tweaks
            gui_slider("_renderer_distance_to_film_plane", &_renderer_distance_to_film_plane, -renderer_n, -renderer_f);
            gui_slider("renderer.persp_distance_to_origin", &renderer.persp_distance_to_origin, 1, 10);
            gui_slider("renderer.theta", &renderer.theta, RAD(-180), RAD(180), true);
            gui_slider("renderer.phi", &renderer.phi, RAD(-90), RAD(90), true);
            gui_slider("renderer.angle_of_view", &renderer.angle_of_view, RAD(1), RAD(178), true);
            gui_slider("renderer_n", &renderer_n, -2, -.1);
            gui_slider("renderer_f", &renderer_f, -20, renderer_n);
            _renderer_distance_to_film_plane = CLAMP(_renderer_distance_to_film_plane, -renderer_n, -renderer_f);
        }

        mat4 C_renderer = camera_get_C(&renderer);
        mat4 P_renderer = _window_get_P_perspective(renderer.angle_of_view, renderer_n, renderer_f, 1); // aspect <- 1
        mat4 V_renderer = inverse(C_renderer);
        rasterize(P_renderer, V_renderer, M, mesh, &color_buffer, &depth_buffer, renderer_n, renderer_f);
        texture_sync_to_GPU(&color_buffer);
        texture_sync_to_GPU(&depth_buffer);

        { // observe
            { // hold tab to check your work
                char *name = "hold to toggle C_observer <- C_renderer";
                static bool clicked;
                bool clicked_this_frame;
                bool released_this_frame = false;
                {
                    clicked_this_frame = gui_button(name, COW_KEY_TAB);
                    if (clicked_this_frame) {
                        clicked = true;
                    }
                    bool selected = (COW1._gui_selected == (void *) name);
                    if (clicked && !selected) {
                        clicked = false;
                        released_this_frame = true;
                    }
                }
                {
                    static Camera3D safe;
                    if (clicked_this_frame) {
                        memcpy(&safe, &observer, sizeof(Camera3D));
                        memcpy(&observer, &renderer, sizeof(Camera3D));
                        _hide_camera_cube = true;
                    }
                    if (released_this_frame) {
                        memcpy(&observer, &safe, sizeof(Camera3D));
                        _hide_camera_cube = false;
                    }
                }
            }

            camera_move(&observer);
            mat4 P_observer = camera_get_P(&observer);
            mat4 V_observer = camera_get_V(&observer);
            mat4 PV_observer = P_observer * V_observer;

            mesh->draw(P_observer, V_observer, M, V3(0, 1, 1));

            { // bespoke widget
                real n = renderer_n;
                real f = renderer_f;
                real hangle = renderer.angle_of_view / 2;
                real r_n = -n * tan(hangle);
                real r_f = -f * tan(hangle);
                { // frustum
                    vec3 frustum_vertex_positions[] = {{r_n,r_n,n},{-r_n,r_n,n},{-r_n,-r_n,n},{r_n,-r_n,n},{r_n,r_n,n},{r_f,r_f,f},{-r_f,r_f,f},{-r_f,-r_f,f},{r_f,-r_f,f},{r_f,r_f,f},{r_f,-r_f,f},{r_n,-r_n,n},{r_f,-r_f,f},{-r_f,-r_f,f},{-r_n,-r_n,n},{-r_f,-r_f,f},{-r_f,r_f,f},{-r_n,r_n,n},{-r_f,r_f,f},{r_f,r_f,f},{r_n,r_n,n},};
                    soup_draw(PV_observer * C_renderer, SOUP_LINE_STRIP, _COUNT_OF(frustum_vertex_positions), frustum_vertex_positions, NULL, monokai.white, 3);
                }
                { // film plane
                    real s = LERP(INVERSE_LERP(_renderer_distance_to_film_plane, -n, -f), r_n, r_f);
                    mesh_draw(
                            P_observer,
                            V_observer,
                            C_renderer * M4_Translation(0, 0, -_renderer_distance_to_film_plane) * M4_Scaling(s),
                            library.meshes.square.num_triangles,
                            library.meshes.square.triangle_indices,
                            library.meshes.square.num_vertices,
                            library.meshes.square.vertex_positions,
                            NULL, NULL, {},
                            library.meshes.square.vertex_texture_coordinates,
                            (!draw_depth_instead_of_color) ? color_buffer.name : depth_buffer.name);
                }

                { // camera
                    library.soups.axes.draw(PV_observer * C_renderer * M4_Scaling(.25));
                    if (!_hide_camera_cube) {
                        library.soups.box.draw(PV_observer * C_renderer * M4_Scaling(.10), .5 * monokai.gray);
                    }
                }
            }
        }
    }

}


int main() {
    APPS {
        APP(hw5a);
    }
    return 0;
}
