// < cow.h >                                                                    
//         \   ^__^                                                             
//          \  (..)\_______           a                                         
//             (oo)\       )\/      quick                                       
//                 ||----w |      and dirty                                     
//                 ||     ||     app library                    james-bern 2022 
//                                                                              
// - cameras use a right-handed coordinate system (point along their negative z)
// - face culling is disabled by default; to enable #define COW_CULL_BACK_FACES 
//   (meshlib is fornow kinda sloppy and will probably break a bit)             


#include "_cow_preamble.cpp"
void xplat_run_to_line() { // debugger entry point
    do_once { xplat_debugbreak(); }
}

#define POINTS GL_POINTS
#define LINES GL_LINES
#define LINE_STRIP GL_LINE_STRIP
#define LINE_LOOP GL_LINE_LOOP
#define TRIANGLES GL_TRIANGLES
#define TRIANGLE_FAN GL_TRIANGLE_FAN
#define TRIANGLE_STRIP GL_TRIANGLE_STRIP
#define QUADS 255
#define QUAD_MESH 254
#define TRIANGLE_MESH 253
#ifndef GL_QUADS
#define GL_QUADS QUADS
#endif
#define XY 2
#define XYZ 3
#define XYZW 4
// #define XY__ 255
#undef RGB // fornow
#define RGB 3
#define RGBA 4

bool initialized;
GLFWwindow *window;
double _macbook_retina_scale; // D:
int widget_active_widget_ID; // fornow


struct {
    #define KEY_TAB GLFW_KEY_TAB
    #define KEY_ESCAPE GLFW_KEY_ESCAPE

    bool key_pressed[512];
    bool key_released[512];
    bool key_held[512];
    bool key_toggle[512];
    bool mouse_left_pressed;
    bool mouse_left_held;
    bool mouse_left_released;
    bool mouse_right_pressed;
    bool mouse_right_held;
    bool mouse_right_released;

    double _mouse_x_Screen;
    double _mouse_y_Screen;
    double _mouse_x_NDC;
    double _mouse_y_NDC;
    double _mouse_dx_Screen;
    double _mouse_dy_Screen;
    double _mouse_dx_NDC;
    double _mouse_dy_NDC;
    double _mouse_wheel_offset;
} input;

struct {
    char *vert = R""""(
        #version 330 core
        layout (location = 0) in vec3 vertex;
        layout (location = 1) in vec4 color;

        out COLOR_POSITION_WRAPPER {
            vec4 color;
            vec4 position;
        } vs_out;

        uniform mat4 transform;
        uniform bool overlay;
        uniform bool has_vertex_colors;
        uniform vec4 fallback_color;

        void main() {
            gl_Position = transform * vec4(vertex, 1);
            vs_out.position = gl_Position;
            vs_out.color = fallback_color;
            if (has_vertex_colors) {
                vs_out.color = color;
            }

            if (overlay) {
                gl_Position.z = -.99 * gl_Position.w; // ?
                vs_out.position = gl_Position;
            }

        }
    )"""";

    char *frag = R""""(
        #version 330 core

        in COLOR_POSITION_WRAPPER {
            vec4 color;
            vec4 position;
        } fs_in;

        out vec4 frag_color;

        void main() {
            fs_in.position; // fornow
            frag_color = fs_in.color;
        }
    )"""";

    char *geom_POINTS = R""""(
        #version 330 core
        layout (points) in;
        layout (triangle_strip, max_vertices = 4) out;
        uniform float aspect;
        uniform float primitive_radius;

        in COLOR_POSITION_WRAPPER {
            vec4 color;
            vec4 position;
        } gs_in[];

        out GS_OUT {
            vec4 color;
            vec4 position;
            vec2 xy;
        } gs_out;

        vec4 position;

        void emit(float x, float y) {
            gs_out.xy = vec2(x, y);
            gl_Position = gs_out.position = (position + sqrt(2) * primitive_radius * vec4(x / aspect, y, 0, 0)) * gl_in[0].gl_Position.w;
            gs_out.color = gs_in[0].color;                                     
            EmitVertex();                                               
        }

        void main() {    
            gs_in[0].position; // fornow
            position = gl_in[0].gl_Position / gl_in[0].gl_Position.w;
            emit(-1, -1);
            emit(1, -1);
            emit(-1, 1);
            emit(1, 1);
            EndPrimitive();
        }  
    )"""";

    char *frag_POINTS = R""""(
        #version 330 core

        in GS_OUT {
            vec4 color;
            vec4 position;
            vec2 xy;
        } fs_in;

        out vec4 frag_color;

        void main() {
            fs_in.color; // fornow
            fs_in.position; // fornow
            fs_in.xy; // fornow
            if (length(fs_in.xy) > 1) {
                discard;
            } else {
                frag_color = fs_in.color;
            }
        }
    )"""";

    char *geom_LINES = R""""(
        #version 330 core
        layout (lines) in;
        layout (triangle_strip, max_vertices = 4) out;
        uniform float aspect;
        uniform float primitive_radius;

        in COLOR_POSITION_WRAPPER {
            vec4 color;
            vec4 position;
        } gs_in[];

        out COLOR_POSITION_WRAPPER {
            vec4 color;
            vec4 position;
        } gs_out;

        void main() {    
            gs_in[0].position; // fornow
            gs_in[1].position; // fornow
            vec4 s = gl_in[0].gl_Position / gl_in[0].gl_Position.w;
            vec4 t = gl_in[1].gl_Position / gl_in[1].gl_Position.w;

            vec4 color_s = gs_in[0].color;
            vec4 color_t = gs_in[1].color;

            vec4 perp = vec4(primitive_radius * vec2(1 / aspect, 1) * normalize(vec2(-1 / aspect, 1) * (t - s).yx), 0, 0);

            gs_out.position = gl_Position = (s + perp) * gl_in[0].gl_Position.w; gs_out.color = color_s; EmitVertex();
            gs_out.position = gl_Position = (s - perp) * gl_in[0].gl_Position.w; gs_out.color = color_s; EmitVertex();
            gs_out.position = gl_Position = (t + perp) * gl_in[1].gl_Position.w; gs_out.color = color_t; EmitVertex();
            gs_out.position = gl_Position = (t - perp) * gl_in[1].gl_Position.w; gs_out.color = color_t; EmitVertex();

            EndPrimitive();
        }  
    )"""";

    int shader_program_POINTS;
    int shader_program_LINES;
    int shader_program_TRIANGLES;

    GLuint VAO[3], VBO[2], EBO[3];
} basic;

struct {
    char *vert = R""""(
        #version 330 core
        layout (location = 0) in vec3 vertex;
        layout (location = 1) in vec3 normal;
        layout (location = 2) in vec3 color;
        layout (location = 3) in vec2 texCoord;

        out BLOCK {
            vec3 position_World;
            vec3 normal_World;
            vec3 color;
            vec2 texCoord;
        } vs_out;

        uniform bool has_vertex_colors;
        uniform vec4 fallback_color;

        uniform bool has_vertex_texCoords;

        uniform mat4 P, V, M, N;

        void main() {
            vec4 tmp = M * vec4(vertex, 1);
            vs_out.position_World = vec3(tmp);
            gl_Position = P * V * tmp;
            vs_out.normal_World = mat3(N) * normal;
            vs_out.color = has_vertex_colors ? color : vec3(fallback_color);
            vs_out.texCoord = has_vertex_texCoords ? texCoord : vec2(0., 0.);
        }
    )"""";

    // world-space blinn-phong, fresnel-schlick with a colorful point light at the eye
    char *frag = R""""(
        #version 330 core

        in BLOCK {
            vec3 position_World;
            vec3 normal_World;
            vec3 color;
            vec2 texCoord;
        } fs_in;

        out vec4 frag_color;

        uniform bool has_vertex_normals;
        uniform bool has_vertex_texCoords;

        uniform sampler2D i_texture;

        uniform vec4 eye_World;


        void main() {
            vec3 world_to_eye = vec3(eye_World) - fs_in.position_World;
            vec3 N = normalize(fs_in.normal_World);
            vec3 E = normalize(world_to_eye);

            vec3 color = fs_in.color;
            float a = 1.;
            if (has_vertex_texCoords) {
                vec4 tmp = texture(i_texture, fs_in.texCoord);
                color = tmp.rgb;
                a = tmp.a;
            }
            if (has_vertex_normals) {
                color *= .8;
                vec3 base = vec3(1);
                float distance = length(world_to_eye);
                float attenuation = 1 / (1 + .02 * distance + .002 * distance * distance);

                vec3 L = E;
                vec3 H = normalize(L + E);
                float F0 = .05;

                float diffuse = max(0, dot(N, L));
                float specular = pow(max(0, dot(N, H)), 100);
                float fresnel = F0 + (1 - F0) * pow(1 - max(0, dot(N, H)), 5);

                color += attenuation * .7 * diffuse * vec3(.7, 1, .7);
                color += attenuation * .7 * specular * vec3(1, .2, .2);
                color += attenuation * .8 * fresnel * vec3(.2, .2, 1);
            }

            frag_color = vec4(color, a);
        }
    )"""";
    int shader_program;
    GLuint VAO, VBO[4], EBO;

    #define FANCY_MAX_NUM_TEXTURES 16
    #define FANCY_MAX_SIZE_TEXTURE_FILENAME 256
    char texture_filenames[FANCY_MAX_NUM_TEXTURES][FANCY_MAX_SIZE_TEXTURE_FILENAME];
    GLuint textures[FANCY_MAX_NUM_TEXTURES];
    int num_textures;
} fancy;

struct {
    double x_curr;
    double y_curr;
    void *selected_widget_ID;
    void *hot_widget_ID;
    double _dx_accumulator;
} imgui;

struct {
    #define MAX_GL_VERTICES 999999
    bool _began;
    bool _called_gl_PV_at_least_once;
    int _num_vertices;
    double _vertex_positions[3 * MAX_GL_VERTICES];
    double _vertex_colors[4 * MAX_GL_VERTICES];
    int _primitive;
    double _size_in_pixels;
    double _PV[16], _M[16];
    bool _multiply_by_M;
    double _color[4];
} gl;



struct Camera2D {
    // this is an "ortho camera" that points at (o_x, o_y)   
    double screen_height_World;
    double o_x;
    double o_y;
};

struct Camera3D {
    // this is an "orbit camera" that points at the origin*  
    //                             *unless _o_x, _o_y nonzero
    //                                                       
    // it is stuck to the surface of a sphere at (theta, phi)
    // these angles can also be interpreted as (yaw, pitch)  
    //                                                       
    // the sphere radius (camera distance) is given by*      
    //     (screen_height_World / 2) / tan(angle_of_view / 2)
    // *we define it this way for compatability with Camera2D
    double distance_to_origin; //sphere_radius
    double angle_of_view;
    double theta; // yaw
    double phi; // pitch
    double _o_x;
    double _o_y;
};

struct BasicMesh { // soup mesh compatible with basic_draw
    int primitive;
    int dimension_of_positions;
    int dimension_of_colors;
    int num_vertices;
    double *vertex_positions;
    double *vertex_colors;
};

struct BasicTriangleMesh3D { // fornow
    int num_vertices;
    vec3 *vertex_positions;
    vec3 *vertex_colors;
};

#ifdef SNAIL_WAS_INCLUDED
struct OrthogonalCoordinateSystem3D {
    //   = [| | | |]   [\ | / |]
    // C = [x y z o] = [- R - t]
    //   = [| | | |]   [/ | \ |]
    //   = [0 0 0 1]   [- 0 - 1]
    mat4 C;
    vec3 x;
    vec3 y;
    vec3 z;
    union {
        vec3 o;
        vec3 t;
    };
    mat4 R; // stored as [\ | / |]
    //           [- R - 0]
    //           [/ | \ |]
    //           [- 0 - 1]
};

struct Texture {
    char *filename;
    int width;
    int height;
    int nrChannels;
    unsigned char *data;
};

struct FancyTriangleMesh3D { // 3D indexed triangle mesh compatible with fancy_draw
    int num_vertices;
    int num_triangles;
    vec3 *vertex_positions;
    vec3 *vertex_normals;
    vec3 *vertex_colors;
    int3 *triangle_indices;
    vec2 *vertex_texCoords;
    char *texture_filename;
};

struct {
    vec3 red    = { 249./255,  38./255, 114./255 };
    vec3 orange = { 253./255, 151./255,  31./255 };
    vec3 yellow = { 255./255, 255./255,  50./255 }; // not the actual monokai yellow cause i don't like it
    vec3 green  = { 166./255, 226./255,  46./255 };
    vec3 blue   = { 102./255, 217./255, 239./255 };
    vec3 purple = { 174./255, 129./255, 255./255 };
    vec3 white  = { 255./255, 255./255, 255./255 };
    vec3 gray   = { 127./255, 127./255, 127./255 };
    vec3 black  = {   0./255,   0./255,   0./255 };
    vec3 brown  = { 123./255,  63./255,   0./255 }; // monokai doesn't actually define a brown
} monokai;

struct {
    int basic_axes_primitive = LINES;
    int basic_axes_vertex_dimension = 3;
    int basic_axes_color_dimension = 3;
    int basic_axes_num_vertices = 6;
    vec3 basic_axes_vertex_positions[6] = { {}, { 1, 0, 0 }, {}, { 0, 1, 0 }, {}, { 0, 0, 1 } };
    vec3 basic_axes_vertex_colors[6] = { monokai.red, monokai.red, monokai.green, monokai.green, monokai.blue, monokai.blue };
    BasicMesh basic_axes = { basic_axes_primitive, basic_axes_vertex_dimension, basic_axes_color_dimension, basic_axes_num_vertices, (double *) basic_axes_vertex_positions, (double *) basic_axes_vertex_colors };

    int basic_tet_primitive = TRIANGLE_MESH;
    int basic_tet_vertex_dimension = 3;
    int basic_tet_color_dimension = 0;
    int basic_tet_num_vertices = 12;
    vec3 basic_tet_vertex_positions[12] = {{0,0,0},{1,0,0},{0,1,0},{0,0,0},{0,1,0},{0,0,1},{0,0,0},{0,0,1},{1,0,0},{1,0,0},{0,1,0},{0,0,1}};
    vec3 *basic_tet_vertex_colors = NULL;
    BasicMesh basic_tet = { basic_tet_primitive, basic_tet_vertex_dimension, basic_tet_color_dimension, basic_tet_num_vertices, (double *) basic_tet_vertex_positions, (double *) basic_tet_vertex_colors };

    int basic_box_primitive = QUAD_MESH;
    int basic_box_vertex_dimension = 3;
    int basic_box_color_dimension = 0;
    int basic_box_num_vertices = 24;
    vec3 basic_box_vertex_positions[24] = {
        {1,1,1},{1,1,-1},{1,-1,-1},{1,-1,1},
        {-1,-1,1}, {-1,-1,-1}, {-1,1,-1}, {-1,1,1},
        {-1,1,1}, {-1,1,-1}, {1,1,-1}, {1,1,1},
        {1,-1,1},{1,-1,-1},{-1,-1,-1},{-1,-1,1},
        {1,1,1},{1,-1,1},{-1,-1,1},{-1,1,1},
        {-1,1,-1}, {-1,-1,-1}, {1,-1,-1}, {1,1,-1},
    };
    vec3 *basic_box_vertex_colors = NULL;
    BasicMesh basic_box = { basic_box_primitive, basic_box_vertex_dimension, basic_box_color_dimension, basic_box_num_vertices, (double *) basic_box_vertex_positions, (double *) basic_box_vertex_colors };

    // just a triangle
    int fancy_tri_num_triangles = 1;
    int3 fancy_tri_triangle_indices[1] = { { 0, 1, 2 } };
    int fancy_tri_num_vertices = 3;
    vec3 fancy_tri_vertex_positions[3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 } };
    vec3 fancy_tri_vertex_normals[3] = { { 0, 0, 1 }, { 0, 0, 1 }, { 0, 0, 1 } };
    FancyTriangleMesh3D fancy_tri = { fancy_tri_num_vertices, fancy_tri_num_triangles, fancy_tri_vertex_positions, fancy_tri_vertex_normals, NULL, fancy_tri_triangle_indices };

    // a square useful for checking if textures are loaded correctly
    int fancy_square_num_triangles = 2;
    int3 fancy_square_triangle_indices[2] = { { 0, 1, 2 }, { 0, 2, 3 } };
    int fancy_square_num_vertices = 4;
    vec3 fancy_square_vertex_positions[4] = { { -1, -1, 0 }, { 1, -1, 0 }, { 1, 1, 0 }, { -1, 1, 0 } };
    vec3 fancy_square_vertex_normals[4] = { { 0, 0, 1 }, { 0, 0, 1 }, { 0, 0, 1 }, { 0, 0, 1 } };
    vec2 fancy_square_texCoords[4] = { { 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 } };
    FancyTriangleMesh3D fancy_square = { fancy_square_num_vertices, fancy_square_num_triangles, fancy_square_vertex_positions, fancy_square_vertex_normals, NULL, fancy_square_triangle_indices, fancy_square_texCoords };
    // [-1, 1]^3;
    // note: we split the box into six regions
    int fancy_box_num_triangles = 12;
    int fancy_box_num_vertices = 24;
    FancyTriangleMesh3D fancy_box = { fancy_box_num_vertices, fancy_box_num_triangles, _meshlib_fancy_box_vertex_positions, _meshlib_fancy_box_vertex_normals, NULL, _meshlib_fancy_box_triangle_indices, _meshlib_fancy_box_vertex_texCoords };

    // unit radius and height; base is centered at origin; y is up
    // note: we split the cone into two regions
    int fancy_cone_num_triangles = 126;
    int fancy_cone_num_vertices = 129;
    FancyTriangleMesh3D fancy_cone = { fancy_cone_num_vertices, fancy_cone_num_triangles, _meshlib_fancy_cone_vertex_positions, _meshlib_fancy_cone_vertex_normals, NULL, _meshlib_fancy_cone_triangle_indices };

    // unit radius and height; base is centered at origin; y is up
    // note: we split the cylinder into three regions
    int fancy_cylinder_num_triangles = 252;
    int fancy_cylinder_num_vertices = 256;
    FancyTriangleMesh3D fancy_cylinder = { fancy_cylinder_num_vertices, fancy_cylinder_num_triangles, _meshlib_fancy_cylinder_vertex_positions, _meshlib_fancy_cylinder_vertex_normals, NULL, _meshlib_fancy_cylinder_triangle_indices };

    // r = 1; centered at origin
    // note: we want smooth normals for the sphere!
    int fancy_sphere_num_triangles = 1280;
    int fancy_sphere_num_vertices = 642;
    FancyTriangleMesh3D fancy_sphere = { fancy_sphere_num_vertices, fancy_sphere_num_triangles, _meshlib_fancy_sphere_vertex_positions, _meshlib_fancy_sphere_vertex_normals, NULL, _meshlib_fancy_sphere_triangle_indices };
} meshlib;
#endif


double util_random_double(double a = 0, double b = 1) {
    return LERP(double(rand()) / RAND_MAX, a, b);
}
int util_random_sign() {
    return util_random_double() < .5 ? -1 : 1;
}
long util_time_in_millis() { // no promises this is even a little bit accurate
    // struct timespec ts;
    // timespec_get(&ts, TIME_UTC);
    // return long(ts.tv_sec * 1000L) + long(ts.tv_nsec / 1000000L);
    using namespace std::chrono;
    milliseconds ms = duration_cast<milliseconds>(system_clock::now().time_since_epoch());
    return (long) ms.count();
}


void linalg_vec3_cross(double *c, double *a, double *b) { // c = a x b
    double tmp[3] = {
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    };
    memcpy(c, tmp, sizeof(tmp));
}
double linalg_vecX_squared_length(int D, double *a) {
    double ret = 0;
    for (int d = 0; d < D; ++d) ret += pow(a[d], 2);
    return ret;
}
double linalg_vecX_squared_distance(int D, double *a, double *b) {
    double ret = 0;
    for (int d = 0; d < D; ++d) ret += pow(a[d] - b[d], 2);
    return ret;
}
void linalg_vecX_normalize(int D, double *a_hat, double *a) { // a_hat = a / |a|
    double L = sqrt(linalg_vecX_squared_length(D, a));
    for (int d = 0; d < D; ++d) a_hat[d] = a[d] / L;
}
#define _4x4(A, i, j) A[4 * (i) + (j)]
void linalg_mat4_times_mat4(double *C, double *A, double *B) { // C = A B
    ASSERT(C);
    ASSERT(A);
    ASSERT(B);

    double tmp[16] = {}; { // allows for e.g. A <- A B;
        for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) for (int k = 0; k < 4; ++k) _4x4(tmp, i, j) += _4x4(A, i, k) * _4x4(B, k, j);
    }
    memcpy(C, tmp, sizeof(tmp));
}
void linalg_mat4_times_vec4_persp_divide(double *b, double *A, double *x) { // b = A x
    ASSERT(b);
    ASSERT(A);
    ASSERT(x);

    double tmp[4] = {}; { // allows for x <- A x
        for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) tmp[i] += _4x4(A, i, j) * x[j];
        if (!IS_ZERO(tmp[3])) {
            ASSERT(!IS_ZERO(x[3]));
            for (int i = 0; i < 4; ++i) tmp[i] /= tmp[3]; 
        }
    }
    memcpy(b, tmp, sizeof(tmp));
}
double linalg_mat4_determinant(double *A) {
    ASSERT(A);
    double A2323 = _4x4(A, 2, 2) * _4x4(A, 3, 3) - _4x4(A, 2, 3) * _4x4(A, 3, 2);
    double A1323 = _4x4(A, 2, 1) * _4x4(A, 3, 3) - _4x4(A, 2, 3) * _4x4(A, 3, 1);
    double A1223 = _4x4(A, 2, 1) * _4x4(A, 3, 2) - _4x4(A, 2, 2) * _4x4(A, 3, 1);
    double A0323 = _4x4(A, 2, 0) * _4x4(A, 3, 3) - _4x4(A, 2, 3) * _4x4(A, 3, 0);
    double A0223 = _4x4(A, 2, 0) * _4x4(A, 3, 2) - _4x4(A, 2, 2) * _4x4(A, 3, 0);
    double A0123 = _4x4(A, 2, 0) * _4x4(A, 3, 1) - _4x4(A, 2, 1) * _4x4(A, 3, 0);
    return _4x4(A, 0, 0) * (_4x4(A, 1, 1) * A2323 - _4x4(A, 1, 2) * A1323 + _4x4(A, 1, 3) * A1223) 
        -  _4x4(A, 0, 1) * (_4x4(A, 1, 0) * A2323 - _4x4(A, 1, 2) * A0323 + _4x4(A, 1, 3) * A0223) 
        +  _4x4(A, 0, 2) * (_4x4(A, 1, 0) * A1323 - _4x4(A, 1, 1) * A0323 + _4x4(A, 1, 3) * A0123) 
        -  _4x4(A, 0, 3) * (_4x4(A, 1, 0) * A1223 - _4x4(A, 1, 1) * A0223 + _4x4(A, 1, 2) * A0123);
}
void linalg_mat4_inverse(double *invA, double *A) {
    ASSERT(invA);
    ASSERT(A);
    double one_over_det = 1 / linalg_mat4_determinant(A);
    double tmp[16] = {}; { // allows for A <- inv(A)
        double A2323 = _4x4(A, 2, 2) * _4x4(A, 3, 3) - _4x4(A, 2, 3) * _4x4(A, 3, 2);
        double A1323 = _4x4(A, 2, 1) * _4x4(A, 3, 3) - _4x4(A, 2, 3) * _4x4(A, 3, 1);
        double A1223 = _4x4(A, 2, 1) * _4x4(A, 3, 2) - _4x4(A, 2, 2) * _4x4(A, 3, 1);
        double A0323 = _4x4(A, 2, 0) * _4x4(A, 3, 3) - _4x4(A, 2, 3) * _4x4(A, 3, 0);
        double A0223 = _4x4(A, 2, 0) * _4x4(A, 3, 2) - _4x4(A, 2, 2) * _4x4(A, 3, 0);
        double A0123 = _4x4(A, 2, 0) * _4x4(A, 3, 1) - _4x4(A, 2, 1) * _4x4(A, 3, 0);
        double A2313 = _4x4(A, 1, 2) * _4x4(A, 3, 3) - _4x4(A, 1, 3) * _4x4(A, 3, 2);
        double A1313 = _4x4(A, 1, 1) * _4x4(A, 3, 3) - _4x4(A, 1, 3) * _4x4(A, 3, 1);
        double A1213 = _4x4(A, 1, 1) * _4x4(A, 3, 2) - _4x4(A, 1, 2) * _4x4(A, 3, 1);
        double A2312 = _4x4(A, 1, 2) * _4x4(A, 2, 3) - _4x4(A, 1, 3) * _4x4(A, 2, 2);
        double A1312 = _4x4(A, 1, 1) * _4x4(A, 2, 3) - _4x4(A, 1, 3) * _4x4(A, 2, 1);
        double A1212 = _4x4(A, 1, 1) * _4x4(A, 2, 2) - _4x4(A, 1, 2) * _4x4(A, 2, 1);
        double A0313 = _4x4(A, 1, 0) * _4x4(A, 3, 3) - _4x4(A, 1, 3) * _4x4(A, 3, 0);
        double A0213 = _4x4(A, 1, 0) * _4x4(A, 3, 2) - _4x4(A, 1, 2) * _4x4(A, 3, 0);
        double A0312 = _4x4(A, 1, 0) * _4x4(A, 2, 3) - _4x4(A, 1, 3) * _4x4(A, 2, 0);
        double A0212 = _4x4(A, 1, 0) * _4x4(A, 2, 2) - _4x4(A, 1, 2) * _4x4(A, 2, 0);
        double A0113 = _4x4(A, 1, 0) * _4x4(A, 3, 1) - _4x4(A, 1, 1) * _4x4(A, 3, 0);
        double A0112 = _4x4(A, 1, 0) * _4x4(A, 2, 1) - _4x4(A, 1, 1) * _4x4(A, 2, 0);

        int i = 0;
        tmp[i++] =  one_over_det * (_4x4(A, 1, 1) * A2323 - _4x4(A, 1, 2) * A1323 + _4x4(A, 1, 3) * A1223);
        tmp[i++] = -one_over_det * (_4x4(A, 0, 1) * A2323 - _4x4(A, 0, 2) * A1323 + _4x4(A, 0, 3) * A1223);
        tmp[i++] =  one_over_det * (_4x4(A, 0, 1) * A2313 - _4x4(A, 0, 2) * A1313 + _4x4(A, 0, 3) * A1213);
        tmp[i++] = -one_over_det * (_4x4(A, 0, 1) * A2312 - _4x4(A, 0, 2) * A1312 + _4x4(A, 0, 3) * A1212);
        tmp[i++] = -one_over_det * (_4x4(A, 1, 0) * A2323 - _4x4(A, 1, 2) * A0323 + _4x4(A, 1, 3) * A0223);
        tmp[i++] =  one_over_det * (_4x4(A, 0, 0) * A2323 - _4x4(A, 0, 2) * A0323 + _4x4(A, 0, 3) * A0223);
        tmp[i++] = -one_over_det * (_4x4(A, 0, 0) * A2313 - _4x4(A, 0, 2) * A0313 + _4x4(A, 0, 3) * A0213);
        tmp[i++] =  one_over_det * (_4x4(A, 0, 0) * A2312 - _4x4(A, 0, 2) * A0312 + _4x4(A, 0, 3) * A0212);
        tmp[i++] =  one_over_det * (_4x4(A, 1, 0) * A1323 - _4x4(A, 1, 1) * A0323 + _4x4(A, 1, 3) * A0123);
        tmp[i++] = -one_over_det * (_4x4(A, 0, 0) * A1323 - _4x4(A, 0, 1) * A0323 + _4x4(A, 0, 3) * A0123);
        tmp[i++] =  one_over_det * (_4x4(A, 0, 0) * A1313 - _4x4(A, 0, 1) * A0313 + _4x4(A, 0, 3) * A0113);
        tmp[i++] = -one_over_det * (_4x4(A, 0, 0) * A1312 - _4x4(A, 0, 1) * A0312 + _4x4(A, 0, 3) * A0112);
        tmp[i++] = -one_over_det * (_4x4(A, 1, 0) * A1223 - _4x4(A, 1, 1) * A0223 + _4x4(A, 1, 2) * A0123);
        tmp[i++] =  one_over_det * (_4x4(A, 0, 0) * A1223 - _4x4(A, 0, 1) * A0223 + _4x4(A, 0, 2) * A0123);
        tmp[i++] = -one_over_det * (_4x4(A, 0, 0) * A1213 - _4x4(A, 0, 1) * A0213 + _4x4(A, 0, 2) * A0113);
        tmp[i++] =  one_over_det * (_4x4(A, 0, 0) * A1212 - _4x4(A, 0, 1) * A0212 + _4x4(A, 0, 2) * A0112);
    }
    memcpy(invA, tmp, sizeof(tmp));
}
void linalg_mat4_transpose(double *AT, double *A) { // AT = A^T
    ASSERT(AT);
    ASSERT(A);
    double tmp[16] = {}; { // allows for A <- transpose(A)
        for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) _4x4(tmp, i, j) = _4x4(A, j, i);
    }
    memcpy(AT, tmp, 16 * sizeof(double));
}



void input_get_mouse_position_and_change_in_position_in_world_coordinates(
        double *PV,
        double *mouse_x_world = NULL,
        double *mouse_y_world = NULL,
        double *mouse_dx_world = NULL,
        double *mouse_dy_world = NULL) {
    ASSERT(PV);
    double World_from_NDC[16] = {};
    linalg_mat4_inverse(World_from_NDC, PV);
    double   xy[4] = { input._mouse_x_NDC,  input._mouse_y_NDC,  0, 1 }; // point
    double dxdy[4] = { input._mouse_dx_NDC, input._mouse_dy_NDC, 0, 0 }; // vector
    linalg_mat4_times_vec4_persp_divide(xy, World_from_NDC, xy);
    linalg_mat4_times_vec4_persp_divide(dxdy, World_from_NDC, dxdy);
    if (mouse_x_world) *mouse_x_world = xy[0];
    if (mouse_y_world) *mouse_y_world = xy[1];
    if (mouse_dx_world) *mouse_dx_world = dxdy[0];
    if (mouse_dy_world) *mouse_dy_world = dxdy[1];
}
#ifdef SNAIL_WAS_INCLUDED
vec2 input_get_mouse_position_in_world_coordinates(mat4 PV) {
    vec2 ret;
    input_get_mouse_position_and_change_in_position_in_world_coordinates(PV.data, &ret.x, &ret.y);
    return ret;
}
vec2 input_get_mouse_change_in_position_in_world_coordinates(mat4 PV) {
    vec2 ret;
    input_get_mouse_position_and_change_in_position_in_world_coordinates(PV.data, NULL, NULL, &ret.x, &ret.y);
    return ret;
}
#endif



void window_set_title(char *title) {
    ASSERT(initialized);
    ASSERT(title);
    glfwSetWindowTitle(window, title);
}
void window_get_dimensions_in_pixels(int *W_in_pixels, int *H_in_pixels) {
    ASSERT(initialized);
    glfwGetFramebufferSize(window, W_in_pixels, H_in_pixels);
}
double window_get_height_in_pixels() {
    ASSERT(initialized);
    int ret, _;
    window_get_dimensions_in_pixels(&_, &ret);
    return ret;
}
double window_get_aspect() {
    int W_in_pixels, H_in_pixels;
    window_get_dimensions_in_pixels(&W_in_pixels, &H_in_pixels);
    return double(W_in_pixels) / H_in_pixels;
}
void window_set_transparency(double a) {
    glfwSetWindowOpacity(window, float(a));
}
void window_get_NDC_from_Screen(double *A) {
    ASSERT(A);

    // [x_NDC] = [1/r_x_Screen            0   0 -1] [x_Screen] = [ x_Screen/r_x_Screen - 1]
    // [y_NDC] = [         0  -1/r_y_Screen   0  1] [y_Screen] = [-y_Screen/r_y_Screen + 1]
    // [0    ] = [         0              0   1  0] [       0] = [                       0]
    // [1    ] = [         0              0   0  1] [       1] = [                       1]

    double r_x_Screen, r_y_Screen; {
        int W_Screen, H_Screen; {
            window_get_dimensions_in_pixels(&W_Screen, &H_Screen);
        }
        r_x_Screen =  W_Screen / 2.;
        r_y_Screen =  H_Screen / 2.;
    }

    memset(A, 0, 16 * sizeof(double));
    _4x4(A, 0, 0) = 1 / r_x_Screen;;
    _4x4(A, 0, 3) = -1;
    _4x4(A, 1, 1) = -1 / r_y_Screen;
    _4x4(A, 1, 3) = 1;
    _4x4(A, 2, 2) = 1;
    _4x4(A, 3, 3) = 1;
}
#ifdef SNAIL_WAS_INCLUDED
vec2 window_get_dimensions_in_pixels() {
    int W, H;
    window_get_dimensions_in_pixels(&W, &H);
    return V2(W, H);
}
#endif



void tform_get_P_perspective(double *P, double angle_of_view, double n = 0, double f = 0, double aspect = 0) {
    if (IS_ZERO(n)) { n = -.1; }
    if (IS_ZERO(f)) { f = -10000; }
    if (IS_ZERO(aspect)) { aspect = window_get_aspect(); }
    ASSERT(P);
    ASSERT(n < 0);
    ASSERT(f < 0);

    // consider a point with coordinates (x, y, -z) in the camera's coordinate system
    //                                                                   where z < 0*
    //                              *recall that the camera's z-axis points backwards

    // 1) imagine projecting the point onto some film plane with height r_y and distance D

    //                r_y                               
    //               -|                                 
    //              - |                                 
    //  angle_y    -  |           y <~ vertex           
    //         \  -   |       -   |                     
    //          |-    |   -       |                     
    //          v     +           |                     
    //         -) -   |           |                     
    //        0-------+-----------+----->               
    //                D          -z                     

    // 2) scale film plane by 1 / r_y to yield NDC film plane (with height 1) and distance Q_y
    // y' is the projected position of vertex y in NDC; i.e., if we can get y', we're done :) 

    //                1 <~ edge of NDC film plane
    //               -|                          
    //              - |                          
    //  angle_y    -  |           y              
    //         \  -   |       -   |              
    //          |-    |   -       |              
    //          v     y'          |              
    //         -) -   |           |              
    //        0-------+-----------+----->        
    //              D / r_y      -z              
    //                ^                          
    //                |                          
    //                cot(angle_y) := Q_y        

    // similar triangles has y' / Q_y = y / -z                     
    //                          => y' = -Q_y * (y / z) (Equation 1)

    // we can repeat this procedure in x      
    // the only difference is Q_x vs. Q_y     
    // -------------------------------------- 
    // cot(angle_x) = D / r_x                 
    // cot(angle_y) = D / r_y                 
    // => r_x cot(angle_x) = r_y cot(angle_y) 
    // recall: aspect := r_x / r_y            
    //  => aspect cot(angle_x) = cot(angle_y) 
    //                  => Q_x = Q_y / aspect.

    // encode Equation 1 (and the variant for x) into a homogeneous matrix equation
    // the third row is a   typical clip plane mapping                             

    // [x'] = [Q_x   0  0  0] [x] = [ Q_x * x] ~> [-Q_x * (x / z)]
    // [y'] = [  0 Q_y  0  0] [y] = [ Q_y * y] ~> [-Q_y * (y / z)]
    // [z'] = [  0   0  a  b] [z] = [  az + b] ~> [      -a - b/z]
    // [ 1] = [  0   0 -1  0] [1] = [      -z] ~> [             1]

    double angle_y = angle_of_view / 2;
    double Q_y = 1 / tan(angle_y);
    double Q_x = Q_y / aspect;

    memset(P, 0, 16 * sizeof(double));
    _4x4(P, 0, 0) = Q_x;
    _4x4(P, 0, 1) = 0; // TERRIBLE PATCH
    _4x4(P, 1, 1) = Q_y;
    _4x4(P, 3, 2) = -1;

    // z'(z) = [-a - b/z]              
    // we want to map [n, f] -> [-1, 1]
    // z'(n) = -a - b/n := -1          
    // z'(f) = -a - b/f :=  1          
    //                                 
    // => a + b/n =  1                 
    //    a + b/f = -1                 
    // => b/n - b/f = 2                
    //                                 
    // => b * (f - n) / (n * f) = 2    
    // => b = (2 * n * f) / (f - n)    
    //                                 
    // => a + (2 * f) / (f - n) = 1    
    // => a = -(n + f) / (f - n)       
    //       = (n + f) / (n - f)       
    _4x4(P, 2, 2) = (n + f) / (n - f);     // a
    _4x4(P, 2, 3) = (2 * n * f) / (f - n); // b
}
void tform_get_P_ortho(double *P, double screen_height_World) {
    ASSERT(P);

    // consider a point with coordinates (x, y, z) in the camera's coordinate system

    // 1) imagine projecting the point onto some film plane with height r_y

    // r_y                                  
    // |                                    
    // |                                    
    // +-----------y                        
    // |           |                        
    // |           |                        
    // +-----------------> minus_z direction

    // 2) scale everything by 1 / r_y to yield NDC film plane (with height 1)

    // 1                                     
    // |                                     
    // |                                     
    // y'----------y / r_y                   
    // |           |                         
    // |           |                         
    // +-----------------> minus_z  direction

    // => y' = y / r_y

    double r_y = screen_height_World / 2;
    double r_x = window_get_aspect() * r_y;

    // [x'] = [1/r_x      0   0  0] [x] = [ x/r_x]
    // [y'] = [    0  1/r_y   0  0] [y] = [ y/r_y]
    // [z'] = [    0      0   a  b] [z] = [az + b]
    // [1 ] = [    0      0   0  1] [1] = [     1]

    memset(P, 0, 16 * sizeof(double));
    _4x4(P, 0, 0) = 1 / r_x;
    _4x4(P, 1, 1) = 1 / r_y;

    // z'(z) = [az + b]                
    // we want to map [n, f] -> [-1, 1]
    // z'(n) = an + b := -1            
    // z'(f) = af + b :=  1            
    //                                 
    // => a * (f - n) = 2              
    //    a = 2 / (f - n)              
    //                                 
    // (2 * f) / (f - n) + b = 1       
    // => b = (n + f) / (n - f)        
    double n = 1000;
    double f = -1000;
    _4x4(P, 2, 2) = 2 / (f - n);
    _4x4(P, 2, 3) = (n + f) / (n - f);

    _4x4(P, 2, 2) *= -1; // sign flip
    _4x4(P, 2, 3) *= -1; //          
    _4x4(P, 3, 3) = 1;
}
void tform_get_PV_ortho(double *PV, double screen_height_World, double eye_x, double eye_y) {
    ASSERT(PV);
    double P[16];
    tform_get_P_ortho(P, screen_height_World);
    double V[16] = { 
        1, 0, 0, -eye_x,
        0, 1, 0, -eye_y,
        0, 0, 1, 0     ,
        0, 0, 0, 1     ,
    };
    linalg_mat4_times_mat4(PV, P, V);
}
void tform_get_PV_hud(double *PV) {
    tform_get_PV_ortho(PV, window_get_height_in_pixels(), 0, 0);
    _4x4(PV, 1, 1) *= -1;
    _4x4(PV, 0, 3) -= 1;
    _4x4(PV, 1, 3) += 1;

    _4x4(PV, 0, 0) *= _macbook_retina_scale;
    _4x4(PV, 1, 1) *= _macbook_retina_scale;
}
double tform_get_distance_to_film_plane(double screen_height_World, double angle_of_view) {
    //                   r_y
    //               -   |  
    //            -      |  
    //         -         |  
    //      - ) theta    |  
    // eye-------------->o  
    //         D            
    double r_y = screen_height_World / 2;
    double theta = angle_of_view / 2;
    return r_y / tan(theta);
}
double tform_get_screen_height_World(double distance_to_film_plane, double angle_of_view) {
    double theta = angle_of_view / 2;
    return 2 * tan(theta) * distance_to_film_plane;
}
#ifdef SNAIL_WAS_INCLUDED
mat4 tform_get_P_perspective(double angle_of_view, double n = 0, double f = 0, double aspect = 0) { mat4 ret; tform_get_P_perspective(ret.data, angle_of_view, n, f, aspect); return ret; }
mat4 tform_get_P_ortho(double screen_height_World) { mat4 ret; tform_get_P_ortho(ret.data, screen_height_World); return ret; }
mat4 tform_get_PV_hud() { mat4 ret; tform_get_PV_hud(ret.data); return ret; }
#endif



void camera_get_PV(Camera2D *camera, double *PV) {
    tform_get_PV_ortho(PV, camera->screen_height_World, camera->o_x, camera->o_y);
}
double camera_get_screen_height_World(Camera3D *camera) {
    return tform_get_screen_height_World(camera->distance_to_origin, camera->angle_of_view);
}
void camera_get_coordinate_system(Camera3D *camera, double *C_out, double *x_out = 0, double *y_out = 0, double *z_out = 0, double *o_out = 0, double *R_out = 0) {
    double camera_o_z = camera->distance_to_origin;

    double C[16]; {
        double T[16] = {
            1, 0, 0, camera->_o_x,
            0, 1, 0, camera->_o_y, // fornow
            0, 0, 1, camera_o_z,
            0, 0, 0, 1,
        };
        double c_x = cos(camera->phi);
        double s_x = sin(camera->phi);
        double c_y = cos(camera->theta);
        double s_y = sin(camera->theta);
        double R_y[16] = {
            c_y , 0, s_y, 0,
            0   , 1, 0  , 0,
            -s_y, 0, c_y, 0,
            0   , 0, 0  , 1,
        };
        double R_x[16] = {
            1,   0,    0, 0, 
            0, c_x, -s_x, 0,
            0, s_x,  c_x, 0,
            0  , 0   , 0, 1,
        };
        linalg_mat4_times_mat4(C,      R_x, T);
        linalg_mat4_times_mat4(C, R_y, C     );
    }

    double xyzo[4][3] = {}; {
        for (int c = 0; c < 4; ++c) for (int r = 0; r < 3; ++r) xyzo[c][r] = _4x4(C, r, c);
    }
    if (C_out) memcpy(C_out, C, sizeof(C));
    {
        if (x_out) memcpy(x_out, xyzo[0], 3 * sizeof(double));
        if (y_out) memcpy(y_out, xyzo[1], 3 * sizeof(double));
        if (z_out) memcpy(z_out, xyzo[2], 3 * sizeof(double));
        if (o_out) memcpy(o_out, xyzo[3], 3 * sizeof(double));
    }
    if (R_out) {
        double tmp[] = {
            C[0], C[1], C[ 2], 0,
            C[4], C[5], C[ 6], 0,
            C[8], C[9], C[10], 0,
            0   , 0   , 0    , 1,
        };
        memcpy(R_out, tmp, sizeof(tmp));
    }
}
void camera_get_P(Camera3D *camera, double *P) {
    tform_get_P_perspective(P, camera->angle_of_view);
    // tform_get_P_ortho(P, camera->screen_height_World);
}
void camera_get_V(Camera3D *camera, double *V) {
    camera_get_coordinate_system(camera, V);
    linalg_mat4_inverse(V, V);
}
void camera_get_PV(Camera3D *camera, double *PV) {
    ASSERT(PV);
    double P[16], V[16];
    camera_get_P(camera, P);
    camera_get_V(camera, V);
    linalg_mat4_times_mat4(PV, P, V);
}
void camera_move(Camera2D *camera, bool disable_pan = false, bool disable_zoom = false) {
    double NDC_from_World[16] = {};
    camera_get_PV(camera, NDC_from_World);
    double x, y, dx, dy;
    input_get_mouse_position_and_change_in_position_in_world_coordinates(NDC_from_World, &x, &y, &dx, &dy);
    if (!disable_pan && input.mouse_right_held) {
        camera->o_x -= dx;
        camera->o_y -= dy;
    }
    else if (!disable_zoom && !IS_ZERO(input._mouse_wheel_offset)) {
        camera->screen_height_World *= (1 - .1 * input._mouse_wheel_offset);

        // zoom while preserving mouse position                
        //                                                     
        // mouse_World' = mouse_World                          
        // mouse_NDC'   = mouse_NDC                            
        //                                                     
        // mouse_NDC'         = mouse_NDC                      
        // P' V' mouse_World'  = mouse_NDC                     
        // P' V' mouse_World   = mouse_NDC                     
        // V' mouse_World     = inv(P') mouse_NDC              
        // mouse_World - eye' = inv(P') mouse_NDC              
        //               eye' = mouse_World - inv(P') mouse_NDC
        //                                    ^----- a -------^

        double mouse_World[4] = { x, y, 0, 1 };
        double a[4] = {};
        double inv_P_prime[16] = {}; {
            Camera2D tmp = { camera->screen_height_World };
            camera_get_PV(&tmp, inv_P_prime);
            linalg_mat4_inverse(inv_P_prime, inv_P_prime);
        }
        double mouse_NDC[4] = {}; {
            linalg_mat4_times_vec4_persp_divide(mouse_NDC, NDC_from_World, mouse_World); 
        }
        linalg_mat4_times_vec4_persp_divide(a, inv_P_prime, mouse_NDC);
        for (int d = 0; d < 2; ++d) (&camera->o_x)[d] = mouse_World[d] - a[d];
    }
}
void camera_move(Camera3D *camera, bool disable_pan = false, bool disable_zoom = false, bool disable_rotate = false) {
    disable_rotate |= (imgui.selected_widget_ID != NULL); // fornow
    disable_rotate |= (widget_active_widget_ID != 0); // fornow
    { // 2D transforms
        Camera2D tmp = { camera_get_screen_height_World(camera), camera->_o_x, camera->_o_y };
        camera_move(&tmp, disable_pan, disable_zoom);
        camera->distance_to_origin = tform_get_distance_to_film_plane(tmp.screen_height_World, camera->angle_of_view);
        camera->_o_x = tmp.o_x;
        camera->_o_y = tmp.o_y;
    }
    if (!disable_rotate && input.mouse_left_held) {
        double fac = 2;
        camera->theta -= fac * input._mouse_dx_NDC;
        camera->phi += fac * input._mouse_dy_NDC;
        camera->phi = CLAMP(camera->phi, -RAD(90), RAD(90));
    }
}
#ifdef SNAIL_WAS_INCLUDED
mat4 camera_get_PV(Camera2D *camera) { mat4 ret; camera_get_PV(camera, ret.data); return ret; }
OrthogonalCoordinateSystem3D camera_get_coordinate_system(Camera3D *camera) {  mat4 C; vec3 x, y, z, o; mat4 R; camera_get_coordinate_system(camera, C.data, x.data, y.data, z.data, o.data, R.data); return { C, x, y, z, o, R }; }
vec3 camera_get_origin(Camera3D *camera) {  vec3 o; camera_get_coordinate_system(camera, NULL, NULL, NULL, NULL, o.data, NULL); return o; }
mat4 camera_get_P(Camera3D *camera) { mat4 ret; camera_get_P(camera, ret.data); return ret; }
mat4 camera_get_V(Camera3D *camera) { mat4 ret; camera_get_V(camera, ret.data); return ret; }
mat4 camera_get_C(Camera3D *camera) { mat4 ret; camera_get_coordinate_system(camera, ret.data); return ret; }
mat4 camera_get_PV(Camera3D *camera) { mat4 ret; camera_get_PV(camera, ret.data); return ret; }
#endif



int shader_compile(char *source, GLenum type) {
    int shader = glCreateShader(type);
    glShaderSource(shader, 1, &source, NULL);
    glCompileShader(shader);
    {
        int success = 0;
        glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
        if (!success) {
            { // log
                char infoLog[512];
                glGetShaderInfoLog(shader, 512, NULL, infoLog);
                printf("%s", infoLog);
            }
            ASSERT(0);
        }
    }
    return shader;
};
int _shader_load_from_file_and_compile(char *filename, GLenum type) {
    FILE *fp = fopen(filename, "r");
    ASSERT(fp);
    fseek(fp, 0L, SEEK_END);
    long lSize = ftell(fp);
    rewind(fp);
    char *buffer = (char *) calloc(1, lSize+1);
    fread(buffer, lSize, 1, fp);

    int ret = shader_compile(buffer, type);

    fclose(fp);
    free(buffer);

    return ret;
}
int shader_build_program(int vertex_shader, int fragment_shader, int geometry_shader = 0) {
    int shader_program = glCreateProgram();
    glAttachShader(shader_program, vertex_shader);
    glAttachShader(shader_program, fragment_shader);
    if (geometry_shader) glAttachShader(shader_program, geometry_shader);
    glLinkProgram(shader_program);
    {
        int success = 0;
        glGetProgramiv(shader_program, GL_LINK_STATUS, &success);
        if (!success) {
            { // log
                char infoLog[512];
                glGetProgramInfoLog(shader_program, 512, NULL, infoLog);
                printf("%s", infoLog);
            }
            ASSERT(0);
        }
    }
    return shader_program;
};
int shader_build_program(char *vertex_shader_source, char *fragment_shader_source, char *geometry_shader_source = 0) {
    int vert = shader_compile(vertex_shader_source, GL_VERTEX_SHADER);
    int frag = shader_compile(fragment_shader_source, GL_FRAGMENT_SHADER);
    int geom = geometry_shader_source ? shader_compile(geometry_shader_source, GL_GEOMETRY_SHADER) : 0;
    return shader_build_program(vert, frag, geom);
}
int shader_get_uniform_location(int ID, char *name) {
    ASSERT(ID);
    int location = glGetUniformLocation(ID, name);
    return location;
}
void shader_set_uniform_double(int ID, char *name, double value) {
    glUniform1f(shader_get_uniform_location(ID, name), (float) value);
}
void shader_set_uniform_int(int ID, char *name, int value) {
    glUniform1i(shader_get_uniform_location(ID, name), value);
}
void shader_set_uniform_bool(int ID, char *name, bool value) {
    glUniform1ui(shader_get_uniform_location(ID, name), value);
}
void shader_set_uniform_vec2(int ID, char *name, double *value) {
    ASSERT(value);
    glUniform2f(shader_get_uniform_location(ID, name), (float) value[0], (float) value[1]);
}
void shader_set_uniform_vec3(int ID, char *name, double *value) {
    ASSERT(value);
    glUniform3f(shader_get_uniform_location(ID, name), (float) value[0], (float) value[1], (float) value[2]);
}
void shader_set_uniform_vec4(int ID, char *name, double *value) {
    ASSERT(value);
    glUniform4f(shader_get_uniform_location(ID, name), (float) value[0], (float) value[1], (float) value[2], (float) value[3]);
}
void shader_set_uniform_mat4(int ID, char *name, double *value) {
    ASSERT(value);
    float as_floats[16] = {}; {
        for (int k = 0; k < 16; ++k) as_floats[k] = float(value[k]);
    }
    glUniformMatrix4fv(shader_get_uniform_location(ID, name), 1, GL_TRUE, as_floats);
}
#ifdef SNAIL_WAS_INCLUDED
void shader_set_uniform_vec2(int ID, char *name, vec2 value) { shader_set_uniform_vec2(ID, name, value.data); }
void shader_set_uniform_vec3(int ID, char *name, vec3 value) { shader_set_uniform_vec3(ID, name, value.data); }
void shader_set_uniform_vec4(int ID, char *name, vec4 value) { shader_set_uniform_vec4(ID, name, value.data); }
void shader_set_uniform_mat4(int ID, char *name, mat4 value) { shader_set_uniform_mat4(ID, name, value.data); }

// fornow
void shader_set_uniform_array_vec3(int ID, char *name, int count, vec3 *value) {
    ASSERT(value);
    float *tmp = (float *) malloc(count * 3 * sizeof(float)); 
    for (int i = 0; i < count; ++i) {
        for (int d = 0; d < 3; ++d) { 
            tmp[3 * i + d] = (float) value[i][d];
        }
    }
    glUniform3fv(shader_get_uniform_location(ID, name), count, tmp);
    free(tmp);
}
#endif



void basic_draw(
        int primitive,
        double *transform,
        int dimension_of_positions,
        int dimension_of_colors,
        int num_vertices,
        double *vertex_positions,
        double *vertex_colors = NULL,
        double r_fallback = 1,
        double g_fallback = 1,
        double b_fallback = 1,
        double a_fallback = 1,
        double size_in_pixels = 0,
        bool overlay = false,
        double r_wireframe = 1,
        double g_wireframe = 1,
        double b_wireframe = 1,
        double a_wireframe = 1) {
    if (num_vertices == 0) { return; } // UPDATE: num_vertices zero is now valid input
    ASSERT(transform);
    ASSERT(dimension_of_positions >= 1 && dimension_of_positions <= 4);
    if (vertex_colors) ASSERT(dimension_of_colors == 3 || dimension_of_colors == 4);
    ASSERT(dimension_of_colors >= 0);
    ASSERT(vertex_positions);

    int mesh_special_case = 0; {
        if (primitive == TRIANGLE_MESH || primitive == QUAD_MESH) {
            basic_draw(primitive == TRIANGLE_MESH ? TRIANGLES : QUADS, transform, dimension_of_positions, dimension_of_colors, num_vertices, vertex_positions, vertex_colors, r_fallback, g_fallback, b_fallback, a_fallback, size_in_pixels, overlay);

            if (primitive == TRIANGLE_MESH) {
                mesh_special_case = 1;
            } else {
                mesh_special_case = 2;
            }

            primitive = LINES;
            vertex_colors = NULL;
            r_fallback = r_wireframe;
            g_fallback = g_wireframe;
            b_fallback = b_wireframe;
            a_fallback = a_wireframe;
        }
    }


    if (IS_ZERO(size_in_pixels)) {
        size_in_pixels = (primitive == POINTS) ? 10 : 5;
    }
    size_in_pixels *= _macbook_retina_scale;

    double fallback_color[4] = { r_fallback, g_fallback, b_fallback, a_fallback };

    glBindVertexArray(basic.VAO[mesh_special_case]);
    int i_attrib = 0;
    auto guarded_push = [&](int buffer_size, void *array, int dim) {
        glDisableVertexAttribArray(i_attrib); // fornow
        if (array) {
            glBindBuffer(GL_ARRAY_BUFFER, basic.VBO[i_attrib]);
            glBufferData(GL_ARRAY_BUFFER, buffer_size, array, GL_DYNAMIC_DRAW);
            glVertexAttribPointer(i_attrib, dim, GL_DOUBLE, 0, 0, NULL);
            glEnableVertexAttribArray(i_attrib);
        }
        ++i_attrib;
    };
    int vvv_size = int(num_vertices * dimension_of_positions * sizeof(double));
    int ccc_size = int(num_vertices * dimension_of_colors * sizeof(double));
    guarded_push(vvv_size, vertex_positions, dimension_of_positions);
    guarded_push(ccc_size, vertex_colors, dimension_of_colors);

    int shader_program = 0; {
        if (primitive == POINTS) {
            shader_program = basic.shader_program_POINTS;
        } else if (primitive == LINES || primitive == LINE_STRIP || primitive == LINE_LOOP) {
            shader_program = basic.shader_program_LINES;
        } else { // including QUADS
            shader_program = basic.shader_program_TRIANGLES;
        }
    }
    ASSERT(shader_program);
    glUseProgram(shader_program);

    shader_set_uniform_double(shader_program, "aspect", window_get_aspect());
    shader_set_uniform_mat4(shader_program, "transform", transform);
    shader_set_uniform_double(shader_program, "primitive_radius", .5 * size_in_pixels / window_get_height_in_pixels());
    shader_set_uniform_bool(shader_program, "has_vertex_colors", vertex_colors != NULL);
    shader_set_uniform_bool(shader_program, "overlay", overlay);
    shader_set_uniform_vec4(shader_program, "fallback_color", fallback_color);

    // glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    if (primitive != QUADS && !mesh_special_case) {
        glDrawArrays(primitive, 0, num_vertices);
    } else {
        // we upload three EBO's _once_               
        // and bind the appropriate one before drawing

        ASSERT(primitive == QUADS || mesh_special_case != 0);

        const int MAX_VERTICES = 1000000;
        ASSERT(num_vertices <= MAX_VERTICES);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, basic.EBO[mesh_special_case]);

        if (primitive == QUADS) {
            primitive = TRIANGLES;
            num_vertices = (num_vertices / 4) * 6;
            static GLuint *indices;
            if (!indices) {
                indices = (GLuint *) malloc(MAX_VERTICES / 4 * 6 * sizeof(GLuint));
                int k = 0;
                for (int i = 0; i < MAX_VERTICES / 4; ++i) {
                    indices[k++] = 4 * i + 2;
                    indices[k++] = 4 * i + 1;
                    indices[k++] = 4 * i + 0;
                    indices[k++] = 4 * i + 3;
                    indices[k++] = 4 * i + 2;
                    indices[k++] = 4 * i + 0;
                }
                glBufferData(GL_ELEMENT_ARRAY_BUFFER, MAX_VERTICES / 4 * 6 * sizeof(GLuint), indices, GL_STATIC_DRAW);
            }
        } else {
            if (mesh_special_case == 1) {
                num_vertices = (num_vertices / 3) * 6;
                static GLuint *indices;
                if (!indices) {
                    indices = (GLuint *) malloc(MAX_VERTICES / 3 * 6 * sizeof(GLuint));
                    int k = 0;
                    for (int i = 0; i < MAX_VERTICES / 3; ++i) {
                        indices[k++] = 3 * i + 0;
                        indices[k++] = 3 * i + 1;
                        indices[k++] = 3 * i + 1;
                        indices[k++] = 3 * i + 2;
                        indices[k++] = 3 * i + 2;
                        indices[k++] = 3 * i + 0;
                    }
                    glBufferData(GL_ELEMENT_ARRAY_BUFFER, MAX_VERTICES / 3 * 6 * sizeof(GLuint), indices, GL_STATIC_DRAW);
                }
            } else {
                ASSERT(mesh_special_case == 2);
                num_vertices = (num_vertices / 4) * 8;
                static GLuint *indices;
                if (!indices) {
                    indices = (GLuint *) malloc(MAX_VERTICES / 4 * 8 * sizeof(GLuint));
                    int k = 0;
                    for (int i = 0; i < MAX_VERTICES / 4; ++i) {
                        indices[k++] = 4 * i + 0;
                        indices[k++] = 4 * i + 1;
                        indices[k++] = 4 * i + 1;
                        indices[k++] = 4 * i + 2;
                        indices[k++] = 4 * i + 2;
                        indices[k++] = 4 * i + 3;
                        indices[k++] = 4 * i + 3;
                        indices[k++] = 4 * i + 0;
                    }
                    glBufferData(GL_ELEMENT_ARRAY_BUFFER, MAX_VERTICES / 4 * 8 * sizeof(GLuint), indices, GL_STATIC_DRAW);
                }
            }
        }
        glDrawElements(primitive, num_vertices, GL_UNSIGNED_INT, NULL);
    }

    // glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}
void basic_draw(
        double *transform,
        BasicMesh mesh,
        double r_fallback,
        double g_fallback,
        double b_fallback,
        double a_fallback,
        double size_in_pixels = 0,
        bool overlay = false,
        double r_wireframe = 1,
        double g_wireframe = 1,
        double b_wireframe = 1,
        double a_wireframe = 1) {
    basic_draw(
            mesh.primitive,
            transform,
            mesh.dimension_of_positions,
            mesh.dimension_of_colors,
            mesh.num_vertices,
            mesh.vertex_positions,
            mesh.vertex_colors ? mesh.vertex_colors : 0,
            r_fallback,
            g_fallback,
            b_fallback,
            a_fallback,
            size_in_pixels,
            overlay,
            r_wireframe,
            g_wireframe,
            b_wireframe,
            a_wireframe);
}
#ifdef SNAIL_WAS_INCLUDED
template<int D_pos, int D_col> void basic_draw(
        int primitive,
        mat4 transform,
        int num_vertices,
        SnailVec<D_pos> *vertex_positions,
        SnailVec<D_col> *vertex_colors,
        double size_in_pixels = 0,
        bool overlay = false,
        vec3 wireframe_color = V3(1, 1, 1)) {
    STATIC_ASSERT(D_pos == 2 || D_pos == 3 || D_pos == 4);
    STATIC_ASSERT(D_col == 3 || D_col == 4);
    basic_draw(
            primitive,
            transform.data,
            D_pos, D_col,
            num_vertices,
            (double *) vertex_positions,
            (double *) vertex_colors,
            0,
            0,
            0,
            0,
            size_in_pixels,
            overlay,
            wireframe_color.r,
            wireframe_color.g,
            wireframe_color.b,
            1);
}
template<int D_pos, int D_col = 3> void basic_draw(
        int primitive,
        mat4 transform,
        int num_vertices,
        SnailVec<D_pos> *vertex_positions,
        SnailVec<D_col> fallback_color = V3(1, 1, 1),
        double size_in_pixels = 0,
        bool overlay = false,
        vec3 wireframe_color = V3(1, 1, 1)) {
    STATIC_ASSERT(D_pos == 2 || D_pos == 3 || D_pos == 4);
    STATIC_ASSERT(D_col == 3 || D_col == 4);
    basic_draw(
            primitive,
            transform.data,
            D_pos,
            D_col,
            num_vertices,
            (double *) vertex_positions,
            NULL,
            fallback_color.r,
            fallback_color.g,
            fallback_color.b,
            D_col == 4 ? fallback_color[3] : 1,
            size_in_pixels,
            overlay,
            wireframe_color.r,
            wireframe_color.g,
            wireframe_color.b,
            1);
}
template<int D_col = 3, int D_col2 = 3> void basic_draw(
        mat4 transform,
        BasicMesh mesh,
        SnailVec<D_col> fallback_color = V3(1, 1, 1),
        double size_in_pixels = 0,
        bool overlay = false,
        SnailVec<3> wireframe_color = V3(1, 1, 1)
        ) {
    STATIC_ASSERT(D_col == 3 || D_col == 4);
    STATIC_ASSERT(D_col2 == 3 || D_col2 == 4);
    basic_draw(transform.data, mesh, fallback_color.r, fallback_color.g, fallback_color.b, D_col == 4 ? fallback_color[3] : 1, size_in_pixels, overlay, wireframe_color.r, wireframe_color.g, wireframe_color.b, D_col2 == 4 ? wireframe_color[3] : 1);
}
void basic_draw(
        int primitive,
        mat4 transform,
        BasicTriangleMesh3D _mesh,
        vec3 fallback_color = V3(1, 1, 1),
        double size_in_pixels = 0,
        bool overlay = false,
        vec3 wireframe_color = V3(1, 1, 1)
        ) {
    ASSERT(primitive == TRIANGLES || primitive == TRIANGLE_MESH);
    BasicMesh mesh = { primitive, 3, 3, _mesh.num_vertices, (double *) _mesh.vertex_positions, (double *) _mesh.vertex_colors };
    basic_draw(transform, mesh, fallback_color, size_in_pixels, overlay, wireframe_color);
}
#endif
void basic_text(
        double *PV,
        char *text,
        double x_world,
        double y_world,
        double z_world = 0,
        double r = 1,
        double g = 1,
        double b = 1,
        double a = 1,
        double font_size_in_pixels = 0,
        double dx_in_pixels = 0,
        double dy_in_pixels = 0,
        bool overlay = true
        ) {
    int W_in_pixels, H_in_pixels; {
        window_get_dimensions_in_pixels(&W_in_pixels, &H_in_pixels);
    }

    if (!PV) {
        static double PV_hud[16];
        tform_get_PV_hud(PV_hud);
        PV = PV_hud;
    }

    if (IS_ZERO(font_size_in_pixels)) font_size_in_pixels = 24;

    font_size_in_pixels *= _macbook_retina_scale;
    dx_in_pixels *= _macbook_retina_scale;
    dy_in_pixels *= _macbook_retina_scale;

    static char buffer[99999]; // ~500 chars
    int num_quads = stb_easy_font_print(0, 0, text, NULL, buffer, sizeof(buffer));
    int num_vertices = 4 * num_quads;
    static double vertex_positions[99999];

    char *read_head = buffer;
    for (int k = 0; k < num_vertices; ++k) {
        for (int d = 0; d < 2; ++d) {
            vertex_positions[2 * k + d] = ((float *) read_head)[d];
        }
        read_head += (3 * sizeof(float) + 4);
    }

    double s_NDC[4] = {}; {
        double s_world[4] = { x_world, y_world, z_world, 1 };
        linalg_mat4_times_vec4_persp_divide(s_NDC, PV, s_world);
    }

    if (IN_RANGE(s_NDC[2], -1, 1)) {

        double transform[16] = {}; {
            // Translation(s_NDC) * app_NDC_from_Screen() * Translation(ds_Screen + dims / 2) * Scaling(size, size);

            double TS[16] = {
                font_size_in_pixels / 12, 0, 0, dx_in_pixels + W_in_pixels / 2,
                0, font_size_in_pixels / 12, 0, dy_in_pixels + H_in_pixels / 2,
                0, 0, 1, 0,
                0, 0, 0, 1,
            };

            double NDC_from_Screen[16] = {}; {
                window_get_NDC_from_Screen(NDC_from_Screen);
            }

            linalg_mat4_times_mat4(transform, NDC_from_Screen, TS);
            for (int d = 0; d < 3; ++d) transform[4 * d + 3] += s_NDC[d];
        }

        basic_draw(QUADS, transform, XY, RGBA, num_vertices, vertex_positions, NULL, r, g, b, a, 0, overlay);
    }
}
#ifdef SNAIL_WAS_INCLUDED
template<int D_pos, int D_color = 3> void basic_text(
        mat4 PV,
        char *text,
        SnailVec<D_pos> s_world,
        SnailVec<D_color> color = V3(1, 1, 1),
        double font_size_in_pixels = 0,
        vec2 ds_in_pixels = {},
        bool overlay = true) {
    STATIC_ASSERT(D_pos == 2 || D_pos == 3 || D_pos == 4);
    STATIC_ASSERT(D_color == 3 || D_color == 4);
    basic_text(PV.data, text, s_world.x, s_world.y, D_pos == 3 ? s_world[2] : 0, color.r, color.g, color.b, D_color == 4 ? color[3] : 1, font_size_in_pixels, ds_in_pixels.x, ds_in_pixels.y, overlay);
}
template<int D_color = 3> void basic_text(
        char *text,
        vec2 s_world,
        SnailVec<D_color> color = V3(1, 1, 1),
        double font_size_in_pixels = 0,
        vec2 ds_in_pixels = {},
        bool overlay = true
        ) {
    STATIC_ASSERT(D_color == 3 || D_color == 4);
    basic_text(NULL, text, s_world.x, s_world.y, 0, color.r, color.g, color.b, D_color == 4 ? color[3] : 1, font_size_in_pixels, ds_in_pixels.x, ds_in_pixels.y, overlay);
}
#endif



bool fancy_texture_find(int *i_texture, char *texture_filename) {
    // FORNOW O(n); TODO hash table
    bool already_loaded = false;
    for ((*i_texture) = 0; (*i_texture) < fancy.num_textures; (*i_texture)++) {
        if (strcmp(texture_filename, fancy.texture_filenames[*i_texture]) == 0) {
            already_loaded = true;
            break;
        }
    }
    return already_loaded;
}
int _fancy_texture_create_load(char *texture_filename) {
    int i_texture;
    ASSERT(fancy.num_textures < FANCY_MAX_NUM_TEXTURES);
    ASSERT(strlen(texture_filename) < FANCY_MAX_SIZE_TEXTURE_FILENAME);
    i_texture = fancy.num_textures++;
    glGenTextures(1, fancy.textures + i_texture);
    strcpy(fancy.texture_filenames[i_texture], texture_filename);

    glActiveTexture(GL_TEXTURE0 + i_texture);
    glBindTexture(GL_TEXTURE_2D, fancy.textures[i_texture]);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);   
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    return i_texture;
}
int _fancy_texture_get_format(int nrChannels) {
    ASSERT(nrChannels == 1 || nrChannels == 3 || nrChannels == 4);
    return (nrChannels == 1) ? GL_RED : (nrChannels == 3) ? GL_RGB : GL_RGBA;
}
int fancy_texture_create(char *texture_filename, int width, int height, int nrChannels, unsigned char *data) {
    int i_texture = _fancy_texture_create_load(texture_filename);
    ASSERT(width > 0);
    ASSERT(height > 0);
    ASSERT(data);
    int format = _fancy_texture_get_format(nrChannels);
    glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, format, GL_UNSIGNED_BYTE, data);
    glGenerateMipmap(GL_TEXTURE_2D);
    return i_texture;
}
void fancy_texture_update(char *texture_filename, int width, int height, int nrChannels, unsigned char *data) {
    ASSERT(width > 0);
    ASSERT(height > 0);
    ASSERT(data);
    int format = _fancy_texture_get_format(nrChannels);

    int i_texture;
    ASSERT(fancy_texture_find(&i_texture, texture_filename));
    glActiveTexture(GL_TEXTURE0 + i_texture);
    glBindTexture(GL_TEXTURE_2D, fancy.textures[i_texture]);

    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, format, GL_UNSIGNED_BYTE, data);
}
int fancy_texture_load(char *texture_filename) {
    int i_texture = _fancy_texture_create_load(texture_filename);

    int width, height, nrChannels;
    unsigned char *data = stbi_load(texture_filename, &width, &height, &nrChannels, 0);
    ASSERT(data);
    ASSERT(nrChannels == 3 || nrChannels == 4);
    int format = (nrChannels == 3) ? GL_RGB : GL_RGBA;
    glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, format, GL_UNSIGNED_BYTE, data);
    glGenerateMipmap(GL_TEXTURE_2D);
    stbi_image_free(data);
    return i_texture;
}
int fancy_texture_create(Texture *texture) {
    return fancy_texture_create(texture->filename, texture->width, texture->height, texture->nrChannels, texture->data);
}
void fancy_texture_update(Texture *texture) {
    fancy_texture_update(texture->filename, texture->width, texture->height, texture->nrChannels, texture->data);
}

// // fancy (circa 1975) draw
// - assumes TRIANGLES, XYZ vertex_positions, RGB vertex_colors
void fancy_draw(
        double *P,
        double *V,
        double *M,
        int num_triangles,
        int *triangle_indices,
        int num_vertices,
        double *vertex_positions,
        double *vertex_normals = NULL,
        double *vertex_colors = NULL,
        double r_fallback = 1,
        double g_fallback = 1,
        double b_fallback = 1,
        double *vertex_texCoords = NULL,
        char *texture_filename = NULL
        ) {
    if (num_triangles == 0) { return; } // NOTE: num_triangles zero is now valid input
    ASSERT(P);
    ASSERT(V);
    ASSERT(M);
    ASSERT(vertex_positions);
    double fallback_color[4] = { r_fallback, g_fallback, b_fallback, 1 };

    int i_texture = -1;
    if (texture_filename) {
        if (!fancy_texture_find(&i_texture, texture_filename)) { // FORNOW load textures on demand
            i_texture = fancy_texture_load(texture_filename);
        }
    }

    glBindVertexArray(fancy.VAO);
    int i_attrib = 0;
    auto guarded_push = [&](int buffer_size, void *array, int dim) {
        glDisableVertexAttribArray(i_attrib); // fornow
        if (array) {
            glBindBuffer(GL_ARRAY_BUFFER, fancy.VBO[i_attrib]);
            glBufferData(GL_ARRAY_BUFFER, buffer_size, array, GL_DYNAMIC_DRAW);
            glVertexAttribPointer(i_attrib, dim, GL_DOUBLE, 0, 0, NULL);
            glEnableVertexAttribArray(i_attrib);
        }
        ++i_attrib;
    };

    int vvv_size = int(num_vertices * 3 * sizeof(double));
    int nnn_size = int(num_vertices * 3 * sizeof(double));
    int ccc_size = int(num_vertices * 3 * sizeof(double));
    int tt_size = int(num_vertices * 2 * sizeof(double));
    guarded_push(vvv_size, vertex_positions, 3);
    guarded_push(nnn_size, vertex_normals, 3);
    guarded_push(ccc_size, vertex_colors, 3);
    guarded_push(tt_size, vertex_texCoords, 2);

    ASSERT(fancy.shader_program);
    glUseProgram(fancy.shader_program);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, fancy.EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3 * num_triangles * sizeof(GLuint), triangle_indices, GL_DYNAMIC_DRAW);

    shader_set_uniform_mat4(fancy.shader_program, "P", P);
    shader_set_uniform_mat4(fancy.shader_program, "V", V);
    shader_set_uniform_mat4(fancy.shader_program, "M", M);
    { // fornow scavenge the camera position from V
        double C[16];
        linalg_mat4_inverse(C, V);
        double eye_World[4] = { C[3], C[7], C[11], 1 };
        shader_set_uniform_vec4(fancy.shader_program, "eye_World", eye_World);
    }
    double N[16]; {
        memcpy(N, M, sizeof(N));
        _4x4(N, 0, 3) = 0;
        _4x4(N, 1, 3) = 0;
        _4x4(N, 2, 3) = 0;
        linalg_mat4_inverse(N, N);
        linalg_mat4_transpose(N, N);
    }
    shader_set_uniform_mat4(fancy.shader_program, "N", N);
    shader_set_uniform_bool(fancy.shader_program, "has_vertex_colors", vertex_colors != NULL);
    shader_set_uniform_bool(fancy.shader_program, "has_vertex_normals", vertex_normals != NULL);
    shader_set_uniform_bool(fancy.shader_program, "has_vertex_texCoords", (vertex_texCoords != NULL) && (texture_filename != NULL));
    shader_set_uniform_int(fancy.shader_program, "i_texture", MAX(0, i_texture));
    shader_set_uniform_vec4(fancy.shader_program, "fallback_color", fallback_color);

    glDrawElements(GL_TRIANGLES, 3 * num_triangles, GL_UNSIGNED_INT, NULL);
}
#ifdef SNAIL_WAS_INCLUDED
void fancy_draw(
        mat4 P,
        mat4 V,
        mat4 M,
        int num_triangles,
        int3 *triangle_indices,
        int num_vertices,
        vec3 *vertex_positions,
        vec3 *vertex_normals = NULL,
        vec3 *vertex_colors = NULL,
        vec3 fallback_color = V3(1, 1, 1),
        vec2 *vertex_texCoords = NULL,
        char *texture_filename = NULL) {
    fancy_draw(P.data, V.data, M.data, num_triangles, (int *) triangle_indices, num_vertices, (double *) vertex_positions, (double *) vertex_normals, (double *) vertex_colors, fallback_color.r, fallback_color.g, fallback_color.b, (double *) vertex_texCoords, texture_filename);
}
void fancy_draw(mat4 P, mat4 V, mat4 M, FancyTriangleMesh3D mesh, vec3 fallback_color = V3(1, 1, 1)) {
    fancy_draw(
            P,
            V,
            M,
            mesh.num_triangles,
            mesh.triangle_indices,
            mesh.num_vertices,
            mesh.vertex_positions,
            mesh.vertex_normals,
            mesh.vertex_colors,
            fallback_color,
            mesh.vertex_texCoords,
            mesh.texture_filename
            );
}
#endif

// 0 if no widget functions active
// otherwise the ID of the active widget function
#define WIDGET_ID_IMGUI 1
#define WIDGET_ID_ANNOTATE 2
#define WIDGET_ID_DRAG 3
#define WIDGET_ID_TRANSLATE 4

void widget_drag(double *PV, int num_vertices, double *vertex_positions, double size_in_pixels = 0, double r = 1, double g = 1, double b = 1, double a = 1) {
    if (widget_active_widget_ID != 0 && widget_active_widget_ID != WIDGET_ID_DRAG) return;
    static double *selected;
    if (selected) { // fornow: allows multiple calls to this function between begin_frame
        bool found = false;
        for (int i = 0; i < num_vertices; ++i) {
            if (selected == vertex_positions + 2 * i) {
                found = true;
                break;
            }
        }
        if (!found) return;
    }
    double *hot = selected;
    if (!selected) {
        {
            for (int i = 0; i < num_vertices; ++i) {
                double *ptr = vertex_positions + 2 * i;
                double tmp[4] = { ptr[0], ptr[1], 0, 1 };
                linalg_mat4_times_vec4_persp_divide(tmp, PV, tmp);
                for (int d = 0; d < 2; ++d) tmp[d] -= (&input._mouse_x_NDC)[d];
                if (linalg_vecX_squared_length(2, tmp) < pow(.02, 2)) { // todo comparison in pixels
                    hot = ptr;
                }
            }
        }
        if (hot) {
            if (input.mouse_left_pressed) {
                selected = hot;
            }
        }
    }

    if (hot || selected) {
        basic_draw(POINTS, PV, XY, RGB, 1, hot, NULL, r, g, b, a, size_in_pixels, true);
    }

    double mouse_xy_World[2];
    input_get_mouse_position_and_change_in_position_in_world_coordinates(PV, mouse_xy_World, mouse_xy_World + 1, NULL, NULL);
    if (input.mouse_left_held && selected) {
        for (int d = 0; d < 2; ++d) *(selected + d) = mouse_xy_World[d];
    } else if (input.mouse_left_released) {
        selected = NULL;
    }

    widget_active_widget_ID = (hot || selected) ? WIDGET_ID_DRAG : 0;
}
#ifdef SNAIL_WAS_INCLUDED
template<int D_color = 3> void widget_drag(mat4 PV, int num_vertices, vec2 *vertex_positions, double size_in_pixels = 0, SnailVec<D_color> color = V3(1, 1, 1)) {
    STATIC_ASSERT(D_color == 3 || D_color == 4);
    return widget_drag(PV.data, num_vertices, (double *) vertex_positions, size_in_pixels, color.r, color.g, color.b, D_color == 4 ? color[3] : 1);
}
#endif



void gl_begin(int primitive, double size_in_pixels = 0) {
    ASSERT(!gl._began);
    gl._began = true;
    gl._primitive = primitive;
    gl._size_in_pixels = size_in_pixels;
    gl._num_vertices = 0;
}
void gl_end() {
    ASSERT(gl._began);
    ASSERT(gl._called_gl_PV_at_least_once);
    gl._began = false;
    basic_draw(gl._primitive,
            gl._PV,
            XYZ,
            RGBA,
            gl._num_vertices,
            gl._vertex_positions,
            gl._vertex_colors,
            0,
            0,
            0,
            0,
            gl._size_in_pixels);
}
void gl_vertex(double x, double y, double z = 0, double w = 1) {
    ASSERT(gl._began);
    double p[4] = { x, y, z, w };
    if (gl._multiply_by_M) linalg_mat4_times_vec4_persp_divide(p, gl._M, p);
    memcpy(gl._vertex_positions + 3 * gl._num_vertices, p, 3 * sizeof(double));
    memcpy(gl._vertex_colors + 4 * gl._num_vertices, gl._color, 4 * sizeof(double));
    ++gl._num_vertices;
}
void gl_color(double r, double g, double b, double a = 1) {
    gl._color[0] = r;
    gl._color[1] = g;
    gl._color[2] = b;
    gl._color[3] = a;
}
void gl_PV(double *PV) {
    gl._called_gl_PV_at_least_once = true;
    memcpy(gl._PV, PV, 16 * sizeof(double));
}
void gl_M(double *M) {
    double zero_matrix__4x4[16] = {};
    memcpy(gl._M, M, 16 * sizeof(double));
    gl._multiply_by_M = (memcmp(gl._M, zero_matrix__4x4, 16 * sizeof(double)) != 0);
}
#ifdef SNAIL_WAS_INCLUDED
void gl_vertex(vec2 s) { gl_vertex(s.x, s.y); }
void gl_vertex(vec3 s) { gl_vertex(s.x, s.y, s.z); }
void gl_vertex(vec4 s) { gl_vertex(s.x, s.y, s.z, s.w); }
void gl_color(vec3 c, double a = 1) { gl_color(c.r, c.g, c.b, a); }
void gl_PV(mat4 PV) { gl_PV(PV.data); }
void gl_M(mat4 M) { gl_PV(M.data); }
#endif







void imgui_begin_frame() {
    imgui.x_curr = 32;
    imgui.y_curr = 48;
}
void _imgui_printf(const char *format, ...) {
    static char text[256] = {};
    va_list arg;
    va_start(arg, format);
    vsnprintf(text, sizeof(text), format, arg);
    va_end(arg);
    basic_text(NULL, text, imgui.x_curr, imgui.y_curr);
    imgui.y_curr += 28;
}
bool imgui_button(char *t, int shortcut_key = 0) {
    double PV[16] = {};
    tform_get_PV_hud(PV);
    double s_mouse[2];
    input_get_mouse_position_and_change_in_position_in_world_coordinates(PV, s_mouse, s_mouse + 1);

    // fornow
    static char text[256];
    if (shortcut_key) {
        if (shortcut_key != KEY_TAB) {
            snprintf(text, sizeof(text), "%s `%c'", t, shortcut_key);
        } else {
            ASSERT(shortcut_key == KEY_TAB);
            snprintf(text, sizeof(text), "%s TAB", t);
        }
    } else {
        strcpy(text, t);
    }
    double L = (2 * stb_easy_font_width(text) + 16); // fornow
    double H = 24;
    double box[8] = {
        imgui.x_curr    , imgui.y_curr    ,
        imgui.x_curr + L, imgui.y_curr    ,
        imgui.x_curr + L, imgui.y_curr + H,
        imgui.x_curr    , imgui.y_curr + H,
    };

    bool hot = IN_RANGE(s_mouse[0], box[0], box[2]) && IN_RANGE(s_mouse[1], box[1], box[5]);

    if (!imgui.selected_widget_ID && ((hot && input.mouse_left_pressed) || input.key_pressed[shortcut_key])) {
        imgui.selected_widget_ID = t;
    }
    if (imgui.selected_widget_ID == t) {
        if (input.mouse_left_released || input.key_released[shortcut_key]) {
            imgui.selected_widget_ID = NULL;
        }
    }

    double r = (imgui.selected_widget_ID != t) ? 0 : .8;
    if (imgui.selected_widget_ID != t) {
        double nudge = SGN(.5 - r) * .1;
        r += nudge; 
        if (hot || input.key_held[shortcut_key]) r += nudge; 
    }
    basic_draw(QUADS, PV, XY, RGB, 4, box, NULL, r, r, r, 1, 0, true);
    basic_draw(LINE_LOOP, PV, XY, RGB, 4, box, NULL, 1, 1, 1, 1, 4, true);
    imgui.x_curr += 8;
    imgui.y_curr += 4;
    _imgui_printf(text);
    imgui.y_curr += 8;
    imgui.x_curr -= 8;

    return (imgui.selected_widget_ID == t) && (input.mouse_left_pressed || input.key_pressed[shortcut_key]);
}
void imgui_checkbox(char *name, bool *t, int shortcut_key = 0) {
    double PV[16] = {};
    tform_get_PV_hud(PV);
    double s_mouse[2];
    input_get_mouse_position_and_change_in_position_in_world_coordinates(PV, s_mouse, s_mouse + 1);
    double L = 16;
    double box[8] = {
        imgui.x_curr    , imgui.y_curr    ,
        imgui.x_curr + L, imgui.y_curr    ,
        imgui.x_curr + L, imgui.y_curr + L,
        imgui.x_curr    , imgui.y_curr + L,
    };
    bool hot = IN_RANGE(s_mouse[0], box[0], box[2]) && IN_RANGE(s_mouse[1], box[1], box[5]);

    if ((hot && input.mouse_left_pressed) || input.key_pressed[shortcut_key]) {
        *t = !(*t);
        if (!imgui.selected_widget_ID) imgui.selected_widget_ID = t;
    }
    if (imgui.selected_widget_ID == t) {
        if (input.mouse_left_released || input.key_released[shortcut_key]) {
            imgui.selected_widget_ID = NULL;
        }
    }

    double r = (!*t) ? 0 : 1;
    if (imgui.selected_widget_ID != t) {
        double nudge = SGN(.5 - r) * .1;
        r += nudge; 
        if (hot || input.key_held[shortcut_key]) r += nudge; 
    }
    basic_draw(QUADS, PV, XY, RGB, 4, box, NULL, r, r, r, 1, 0, true);
    basic_draw(LINE_LOOP, PV, XY, RGB, 4, box, NULL, 1, 1, 1, 1, 4, true);
    imgui.x_curr += 2 * L;
    if (shortcut_key) {
        if (shortcut_key != KEY_TAB) {
            _imgui_printf("%s `%c'", name, shortcut_key);
        } else {
            _imgui_printf("%s `TAB'", name);
        }
    } else {
        _imgui_printf(name);
    }
    imgui.x_curr -= 2 * L;
}
void imgui_readout(char *name, int *t) {
    if (!name) name = "";
    char *join = (char *)((name) ? " " : "");
    _imgui_printf("%s%s%d", name, join, *t);
}
void imgui_readout(char *name, double *t) {
    if (!name) name = "";
    char *join = (char *)((name) ? " " : "");
    _imgui_printf("%s%s%lf", name, join, *t);
}
void imgui_readout(char *name, Camera2D *t) {
    if (!name) name = "";
    char *join0 = (char *)((name) ? "->" : "");
    char *join1 = " ";
    #define Q(field) _imgui_printf("%s%s%s%s%lf", name, join0, XSTR(field), join1, t->field);
    Q(screen_height_World);
    Q(o_x);
    Q(o_y);
    #undef Q
}
void imgui_readout(char *name, Camera3D *t) {
    if (!name) name = "";
    char *join0 = (char *)((name) ? "->" : "");
    char *join1 = " ";
    #define Q(field) _imgui_printf("%s%s%s%s%lf", name, join0, XSTR(field), join1, t->field);
    #define R(field) _imgui_printf("%s%s%s%s%lf (%d deg)", name, join0, XSTR(field), join1, t->field, int(round(DEG(t->field))));
    Q(distance_to_origin);
    R(angle_of_view);
    R(theta);
    R(phi);
    Q(_o_x);
    Q(_o_y);
    #undef Q
    #undef R
}
void _imgui_slider(char *text, void *t, bool is_int, double *t_copy, double a, double b, bool display_in_degrees = false) {
    imgui.y_curr += 8;
    double PV[16] = {};
    tform_get_PV_hud(PV);
    double s_mouse[2];
    input_get_mouse_position_and_change_in_position_in_world_coordinates(PV, s_mouse, s_mouse + 1);
    double w = 166;
    double band[4] = { imgui.x_curr, imgui.y_curr, imgui.x_curr + w, imgui.y_curr };
    double s_dot[2] = { LERP(INVERSE_LERP(*t_copy, a, b), band[0], band[2]), band[1] };

    if (widget_active_widget_ID == 0 || widget_active_widget_ID == WIDGET_ID_IMGUI) {
        bool is_near = linalg_vecX_squared_distance(2, s_dot, s_mouse) < _macbook_retina_scale * 16;
        if (is_near) {
            imgui.hot_widget_ID = t;
            widget_active_widget_ID = WIDGET_ID_IMGUI;
        }
        if (imgui.hot_widget_ID == t && !is_near) {
            imgui.hot_widget_ID = 0;
            if (imgui.selected_widget_ID != t) widget_active_widget_ID = 0;
        }
    }
    if (!imgui.selected_widget_ID && (imgui.hot_widget_ID == t) && input.mouse_left_pressed) {
        widget_active_widget_ID = WIDGET_ID_IMGUI;
        imgui.selected_widget_ID = t;
    }
    if (imgui.selected_widget_ID == t) {
        if (input.mouse_left_held) {
            *t_copy = LERP(CLAMP(INVERSE_LERP(s_mouse[0], band[0], band[2]), 0, 1), a, b);
        }
        if (input.mouse_left_released) {
            if (imgui.hot_widget_ID != t) widget_active_widget_ID = 0;
            imgui.selected_widget_ID = 0;
        }
    }
    basic_draw(LINES, PV, XY, RGB, 2, band, NULL, .6, .6, .6, 1, 6, true);
    double r = (imgui.selected_widget_ID == t) ? 1 : (imgui.hot_widget_ID == t) ? .9 : .8;
    basic_draw(POINTS, PV, XY, RGB, 1, s_dot, NULL, r, r, r, 1, (imgui.hot_widget_ID == t && imgui.selected_widget_ID != t) ? 12 : 10, true);
    imgui.y_curr -= 8;
    imgui.x_curr += w + 16;
    if (is_int) {
        imgui_readout(text, (int *) t);
    } else {
        if (!display_in_degrees) {
            imgui_readout(text, (double *) t);
        } else {
            int tmp = (int) round(DEG(*((double *) t)));
            _imgui_printf("%s %d deg", text, tmp);
        }
    }
    imgui.x_curr -= w + 16;
}
void imgui_slider(char *name, int *t, int a, int b, int j = 0, int k = 0, bool loop = false) {
    double tmp = double(*t);
    static char text[256]; {
        if (!j && !k) {
            snprintf(text, sizeof(text), "%s", name);
        } else {
            if (k == KEY_TAB) {
                ASSERT(!j);
                snprintf(text, sizeof(text), "%s TAB", name);
            } else {
                snprintf(text, sizeof(text), "%s `%c' `%c'", name, j ? j : '~', k ? k : '~');
            }
        }
    }
    _imgui_slider(text, t, true, &tmp, a, b);
    *t = int(.5 + tmp);
    if (input.key_pressed[k]) ++(*t);
    if (input.key_pressed[j]) --(*t);
    if (input.key_pressed[k] || input.key_pressed[j]) {
        *t = (!loop) ? CLAMP(*t, a, b) : a + MODULO(*t - a, (b + 1) - a);
    }
}
void imgui_slider(char *name, double *t, double a, double b, bool display_in_degrees = false) {
    _imgui_slider(name, t, false, t, a, b, display_in_degrees);
}



#ifdef SNAIL_WAS_INCLUDED
vec3 color_get_kelly(int i) {
    static vec3 _kelly_colors[]={{255./255,179./255,0./255},{128./255,62./255,117./255},{255./255,104./255,0./255},{166./255,189./255,215./255},{193./255,0./255,32./255},{206./255,162./255,98./255},{129./255,112./255,102./255},{0./255,125./255,52./255},{246./255,118./255,142./255},{0./255,83./255,138./255},{255./255,122./255,92./255},{83./255,55./255,122./255},{255./255,142./255,0./255},{179./255,40./255,81./255},{244./255,200./255,0./255},{127./255,24./255,13./255},{147./255,170./255,0./255},{89./255,51./255,21./255},{241./255,58./255,19./255},{35./255,44./255,22./255}};
    return _kelly_colors[MODULO(i, NELEMS(_kelly_colors))];
}
vec3 color_rainbow_swirl(double t) {
    #define Q(o) (.5 + .5 * cos(6.28 * ((o) - t)))
    return { Q(0), Q(.33), Q(-.33) };
    #undef Q
}
#endif



double _callback_scaling_factor() {
    return input.key_held[GLFW_KEY_LEFT_SHIFT] ? .1 : 1;
}
void callback_key(GLFWwindow *, int key, int, int action, int) {
    if (key < 0) { return; }
    if (action == GLFW_PRESS) {
        input.key_pressed[key] = true;
        input.key_held[key] = true;
        input.key_toggle[key] = !input.key_toggle[key];
    } else if (action == GLFW_RELEASE) {
        input.key_released[key] = true;
        input.key_held[key] = false;
    }
}
void callback_cursor_position(GLFWwindow *, double xpos, double ypos) {
    double tmp_mouse_x_NDC = input._mouse_x_NDC;
    double tmp_mouse_y_NDC = input._mouse_y_NDC;
    double tmp_mouse_x_Screen = input._mouse_x_Screen;
    double tmp_mouse_y_Screen = input._mouse_y_Screen;

    { // macbook retina nonsense
        xpos *= _macbook_retina_scale;
        ypos *= _macbook_retina_scale;
    }

    input._mouse_x_Screen = xpos;
    input._mouse_y_Screen = ypos;
    double s_NDC[4] = {}; {
        double NDC_from_Screen[16] = {}; {
            window_get_NDC_from_Screen(NDC_from_Screen);
        }
        double s_Screen[4] = { xpos, ypos, 0, 1 };
        linalg_mat4_times_vec4_persp_divide(s_NDC, NDC_from_Screen, s_Screen);
    }
    input._mouse_x_NDC = s_NDC[0];
    input._mouse_y_NDC = s_NDC[1];
    // callback may be called multiple times per frame
    input._mouse_dx_Screen += _callback_scaling_factor() * (input._mouse_x_Screen - tmp_mouse_x_Screen);
    input._mouse_dy_Screen += _callback_scaling_factor() * (input._mouse_y_Screen - tmp_mouse_y_Screen);
    input._mouse_dx_NDC += _callback_scaling_factor() * (input._mouse_x_NDC - tmp_mouse_x_NDC);
    input._mouse_dy_NDC += _callback_scaling_factor() * (input._mouse_y_NDC - tmp_mouse_y_NDC);
}
void callback_mouse_button(GLFWwindow *, int button, int action, int) {
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        if (action == GLFW_PRESS) { 
            input.mouse_left_pressed = true;
            input.mouse_left_held = true;
        } else if (action == GLFW_RELEASE) { 
            input.mouse_left_released = true;
            input.mouse_left_held = false;
        }
    } else if (button == GLFW_MOUSE_BUTTON_RIGHT) {
        if (action == GLFW_PRESS) { 
            input.mouse_right_pressed = true;
            input.mouse_right_held = true;
        } else if (action == GLFW_RELEASE) { 
            input.mouse_right_held = false;
        }
    }
}
void callback_scroll(GLFWwindow *, double, double yoffset) {
    input._mouse_wheel_offset += _callback_scaling_factor() * yoffset;
}
void callback_framebuffer_size(GLFWwindow *, int width, int height) {
    glViewport(0, 0, width, height);
}



void poll_input() {
    ASSERT(initialized);
    memset(input.key_pressed, 0, sizeof(input.key_pressed));
    memset(input.key_released, 0, sizeof(input.key_released));
    input.mouse_left_pressed = false;
    input.mouse_left_released = false;
    input.mouse_right_pressed = false;
    input._mouse_dx_Screen = 0;
    input._mouse_dy_Screen = 0;
    input._mouse_dx_NDC = 0;
    input._mouse_dy_NDC = 0;
    input._mouse_wheel_offset = 0;
    glfwPollEvents();
    // make e.g. key_*['j'] the same as key_*['J']
    for (int i = 0; i < 26; ++i) {
        input.key_pressed ['a' + i] = input.key_pressed ['A' + i];
        input.key_held    ['a' + i] = input.key_held    ['A' + i];
        input.key_released['a' + i] = input.key_released['A' + i];
        input.key_toggle  ['a' + i] = input.key_toggle  ['A' + i];
    }
}
void swap_draw_buffers() {
    ASSERT(initialized);
    glfwSwapBuffers(window);
}
void clear_draw_buffer(double r, double g, double b, double a) {
    ASSERT(initialized);
    glClearColor(float(r), float(g), float(b), float(a));
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
}
bool begin_frame(double r = 0, double g = 0, double b = 0, double a = 0) {
    poll_input();
    swap_draw_buffers();
    clear_draw_buffer(r, g, b, a);
    imgui_begin_frame();
    { // framerate
        static int measured_fps;
        // request uncapped framerate 
        if (input.key_pressed['/']) {
            glfwSwapInterval(!input.key_toggle['/']);
        }
        // grab and smooth fps
        {
            const int N_MOVING_WINDOW = 5;
            static long prev_stamps[N_MOVING_WINDOW];
            long stamp = util_time_in_millis();
            measured_fps = (int) round(N_MOVING_WINDOW / (double(stamp - prev_stamps[N_MOVING_WINDOW - 1]) / 1000.));
            for (int i = N_MOVING_WINDOW - 1; i >= 1; --i) {
                prev_stamps[i] = prev_stamps[i - 1];
            }
            prev_stamps[0] = stamp;
        }
        // display fps
        if (input.key_toggle['\\']) {
            static int display_fps;
            static long stamp = util_time_in_millis();
            if (util_time_in_millis() - stamp > 166) {
                stamp = util_time_in_millis();
                display_fps = measured_fps;
                // printf("fps: %d\n", display_fps);
            }
            char text[256] = {};
            snprintf(text, sizeof(text), "fps: %d", display_fps);
            basic_text(NULL, text, 0, 0, 0, (display_fps < 45) ? 1 : 0, (display_fps > 30) ? 1 : 0, 0);
        }
    }
    return !(input.key_pressed['Q'] || glfwWindowShouldClose(window));
}
void enable_floating_point_exceptions() {
    // ! can't trap invalid on apple silicon (glfw throws)
    // ! apple intel untested (i don't have one)          

    #if defined(unix) || defined(__unix__) || defined(__unix) // ubuntu
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    #elif defined(__APPLE__) || defined(__MACH__) // mac

    #if defined(__arm__) || defined(__arm64__) // Apple Silicon
    fenv_t env;
    fegetenv(&env);
    env.__fpcr = env.__fpcr | __fpcr_trap_divbyzero /* | __fpcr_trap_overflow | __fpcr_trap_invalid */;
    fesetenv(&env);
    #else // Intel
    // _mm_setcsr(_MM_MASK_MASK & ~(_MM_MASK_DIV_ZERO | _MM_MASK_INVALID | _MM_MASK_OVERFLOW));
    #endif

    #elif defined(WIN32) || defined(_WIN64) // windows
    _clearfp();
    _controlfp(_controlfp(0, 0) & ~(_EM_ZERODIVIDE | _EM_INVALID | _EM_OVERFLOW), _MCW_EM);
    #endif
}
void init(bool transparent_framebuffer = false, char *window_title = 0, int screen_height_in_pixels = 540, int window_top_left_init_x = 0, int window_top_left_init_y = 30) {

    if (initialized) {
        if (window_title) window_set_title(window_title);
        {
            double _mouse_x_Screen = input._mouse_x_Screen;
            double _mouse_y_Screen = input._mouse_y_Screen;
            double _mouse_x_NDC = input._mouse_x_NDC;
            double _mouse_y_NDC = input._mouse_y_NDC;
            memset(&input, 0, sizeof(input));
            input._mouse_x_Screen = _mouse_x_Screen;
            input._mouse_y_Screen = _mouse_y_Screen;
            input._mouse_x_NDC = _mouse_x_NDC;
            input._mouse_y_NDC = _mouse_y_NDC;
        }
        widget_active_widget_ID = 0;
        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
        fancy.num_textures = 0;
        memset(fancy.textures, 0, sizeof(fancy.textures));
        memset(&gl, 0, sizeof(gl));
        return;
    }

    // floating point exceptions
    #ifdef COW_CRASH_ON_FLOATING_POINT_EXCEPTIONS
    enable_floating_point_exceptions();
    #endif

    setvbuf(stdout, NULL, _IONBF, 0); // don't buffer printf
    srand(0);
    // srand((unsigned int) time(NULL));

    { // glfw, gl
        ASSERT(glfwInit());

        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
        if (transparent_framebuffer) glfwWindowHint(GLFW_TRANSPARENT_FRAMEBUFFER, GLFW_TRUE);

        // glfwWindowHint(GLFW_SAMPLES, 4); // multisampling

        window = glfwCreateWindow(16 * screen_height_in_pixels / 9, screen_height_in_pixels, (window_title) ? window_title : "cow.cpp! :D", NULL, NULL);
        if (!window) {
            printf("[cow] something's gone wonky; if you weren't just messing with init(...) or something, please try restarting your computer and try again.\n");
            ASSERT(0);
        }

        glfwSetWindowPos(window, window_top_left_init_x, window_top_left_init_y);
        glfwMakeContextCurrent(window);
        glfwSetFramebufferSizeCallback(window, callback_framebuffer_size);
        glfwSetFramebufferSizeCallback(window, callback_framebuffer_size);

        #if defined(WIN32) || defined(_WIN64) // windows ///////////////////////////////
        gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);
        #endif /////////////////////////////////////////////////////////////////////////

        glEnable(GL_DEPTH_TEST);
        glDepthMask(GL_TRUE);
        glDepthFunc(GL_LEQUAL);
        glDepthRange(0.0f, 1.0f);

        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        #ifdef COW_CULL_BACK_FACES
        glEnable(GL_CULL_FACE);  
        #else
        glDisable(GL_CULL_FACE);  
        #endif

        glfwSwapInterval(1);

        glfwSetKeyCallback(window, callback_key);
        glfwSetCursorPosCallback(window, callback_cursor_position);
        glfwSetMouseButtonCallback(window, callback_mouse_button);
        glfwSetScrollCallback(window, callback_scroll);

        initialized = true;
        clear_draw_buffer(0, 0, 0, 0);
    }

    { // _macbook_retina_scale D:
        int num, den, _;
        glfwGetFramebufferSize(window, &num, &_);
        glfwGetWindowSize(window, &den, &_);
        _macbook_retina_scale = num / den;
    }

    { // basic
        basic.shader_program_POINTS = shader_build_program(basic.vert, basic.frag_POINTS, basic.geom_POINTS);
        basic.shader_program_LINES = shader_build_program(basic.vert, basic.frag, basic.geom_LINES);
        basic.shader_program_TRIANGLES = shader_build_program(basic.vert, basic.frag);
        glGenVertexArrays(NELEMS(basic.VAO), basic.VAO);
        glGenBuffers(NELEMS(basic.VBO), basic.VBO);
        glGenBuffers(NELEMS(basic.EBO), basic.EBO);
    }

    { // fancy
        fancy.shader_program = shader_build_program(fancy.vert, fancy.frag);
        glGenVertexArrays(1, &fancy.VAO);
        glGenBuffers(NELEMS(fancy.VBO), fancy.VBO);
        glGenBuffers(1, &fancy.EBO);
        stbi_set_flip_vertically_on_load(true); // stb_image
    }

}






void hello() {
    init();

    Camera2D camera = { 3 };

    bool playing = false;
    vec3 color = { 1, 1, 1 };
    vec2 vertex_positions[3] = {
        { cos(RAD( 90)), sin(RAD( 90)) },
        { cos(RAD(210)), sin(RAD(210)) },
        { cos(RAD(330)), sin(RAD(330)) },
    };

    while (begin_frame()) {
        // xplat_run_to_line();

        camera_move(&camera);
        mat4 PV = camera_get_PV(&camera);

        imgui_checkbox("playing", &playing, 'p');
        if (playing || input.key_pressed['.']) {
            // hmm...
            double fac = .01;
            vertex_positions[0].x -= fac * vertex_positions[0].y;
            vertex_positions[0].y += fac * vertex_positions[0].x;
            vertex_positions[1].x -= fac * vertex_positions[1].y;
            vertex_positions[1].y += fac * vertex_positions[1].x;
            vertex_positions[2].x -= fac * vertex_positions[2].y;
            vertex_positions[2].y += fac * vertex_positions[2].x;
        }

        imgui_slider("r", &color.r, 0, 1);
        imgui_slider("g", &color.g, 0, 1);
        imgui_slider("b", &color.b, 0, 1);
        imgui_readout("x0", &vertex_positions[0].x); imgui_readout("y0", &vertex_positions[0].y);
        imgui_readout("x1", &vertex_positions[1].x); imgui_readout("y1", &vertex_positions[1].y);
        imgui_readout("x2", &vertex_positions[2].x); imgui_readout("y2", &vertex_positions[2].y);
        basic_draw(TRIANGLES, PV, 3, vertex_positions, color);
        widget_drag(PV, 3, vertex_positions);
    }                                                              
}                                                                  

