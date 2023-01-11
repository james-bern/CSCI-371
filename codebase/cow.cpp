#ifndef COW_CPP
#define COW_CPP

////////////////////////////////////////////////////////////////////////////////
// #include "cow.h"/////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#define CUTE_SOUND_IMPLEMENTATION
#include "ext/cute_sound.h"
#include "ext/stb_easy_font.h"
#define STB_IMAGE_IMPLEMENTATION
#include "ext/stb_image.h"
#include "ext/jo_mpeg.h"

#if defined(unix) || defined(__unix__) || defined(__unix)
#define GL_GLEXT_PROTOTYPES
#include <GLFW/glfw3.h>
#elif defined(__APPLE__) || defined(__MACH__)
#define GL_SILENCE_DEPRECATION
#define GLFW_INCLUDE_GL_COREARB
#include <OpenGL/gl3.h>
#include <GLFW/glfw3.h>
#elif defined(WIN32) || defined(_WIN32) || defined(_WIN64)
#include "ext/glad.c"
#include "ext/glfw3.h"
#define GLFW_EXPOSE_NATIVE_WIN32
#define GLFW_EXPOSE_NATIVE_WGL
#define GLFW_NATIVE_INCLUDE_NONE
#include "ext/glfw3native.h"
#else
#pragma message("[cow] operating system not recognized")
#endif

#include <cmath>
#include <time.h>
#include <chrono>

#define COW_MOUSE_OWNER_NONE 0
#define COW_MOUSE_OWNER_GUI 1
#define COW_MOUSE_OWNER_WIDGET 2
#define COW_KEY_TAB GLFW_KEY_TAB
#define COW_KEY_ESCAPE GLFW_KEY_ESCAPE
#define COW_KEY_ARROW_LEFT GLFW_KEY_LEFT
#define COW_KEY_ARROW_RIGHT GLFW_KEY_RIGHT
#undef GLFW_KEY_LEFT_SHIFT  // use cow.key_shift_held instead
#undef GLFW_KEY_RIGHT_SHIFT // use cow.key_shift_held instead
#define SOUND_MAX_DIFFERENT_FILES 32
#define SOUND_MAX_FILENAME_LENGTH 64
#define ITRI_MAX_NUM_TEXTURES 32
#define ITRI_MAX_FILENAME_LENGTH 64

typedef unsigned char u8;
typedef unsigned int u32;

// use these ONLY in cow.cpp in structs where you would really like to store a snail type
#ifdef SNAIL_CPP
typedef vec2 _vec2;
typedef vec3 _vec3;
typedef vec4 _vec4;
typedef mat4 _mat4;
#else
typedef real _vec2[2];
typedef real _vec3[3];
typedef real _vec4[4];
typedef real _mat4[16];
#endif

struct CW_USER_FACING_CONFIG {
    int hotkeys_app_next = COW_KEY_ARROW_RIGHT;
    int hotkeys_app_prev = COW_KEY_ARROW_LEFT;
    int hotkeys_app_menu = '0';
    int hotkeys_gui_hide = '-';

    bool tweaks_soup_draw_with_rounded_corners_for_all_line_primitives = false;
    bool tweaks_ASSERT_crashes_the_program_without_you_having_to_press_Enter = false;
    bool tweaks_record_raw_then_encode_everything_WARNING_USES_A_LOT_OF_DISK_SPACE = false;
    real tweaks_size_in_pixels_soup_draw_defaults_to_if_you_pass_0_for_size_in_pixels = 12.0;

    _vec2 default_window_position = { 0.0, 30.0 };
    _vec2 default_window_size = { 1280.0, 720.0 };
    _vec3 _default_window_color_rgb = { 0.0, 0.0, 0.0 };
    real default_window_color_a = 1.0;
    bool default_window_floating = false;
    bool default_window_decorated = true;
};

struct C2_READONLY_USER_FACING_DATA {
    bool key_pressed[512];
    bool key_held[512];
    bool key_released[512];
    bool key_shift_held;

    bool mouse_left_pressed;
    bool mouse_left_held;
    bool mouse_left_released;
    bool mouse_right_pressed;
    bool mouse_right_held;
    bool mouse_right_released;
    real mouse_wheel_offset;
    _vec2 mouse_position_Screen;
    _vec2 mouse_position_NDC;
    _vec2 mouse_change_in_position_Screen;
    _vec2 mouse_change_in_position_NDC;

    int _mouse_owner;

    _mat4 I = { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0 };
    _mat4 NDC_from_Screen;
};

struct CX_INTERNAL_CONSTANTS {
    char *_soup_vert = R"(
        #version 330 core
        layout (location = 0) in vec3 vertex;
        layout (location = 1) in vec4 color;

        out BLOCK {
            vec4 color;
        } vs_out;

        uniform mat4 PVM;
        uniform bool force_draw_on_top;
        uniform bool has_vertex_colors;
        uniform vec4 color_if_vertex_colors_is_NULL;

        void main() {
            gl_Position = PVM * vec4(vertex, 1);
            vs_out.color = color_if_vertex_colors_is_NULL;

            if (has_vertex_colors) {
                vs_out.color = color;
            }

            if (force_draw_on_top) {
                gl_Position.z = -.99 * gl_Position.w; // ?
            }
        }
    )";

    char *_soup_geom_POINTS = R"(
        #version 330 core
        layout (points) in;
        layout (triangle_strip, max_vertices = 4) out;
        uniform float aspect;
        uniform float primitive_radius_NDC;

        in BLOCK {
            vec4 color;
        } gs_in[];

        out GS_OUT {
            vec4 color;
            vec2 xy;
        } gs_out;

        vec4 _position;

        void emit(float x, float y) {
            gs_out.xy = vec2(x, y);
            gl_Position = (_position + primitive_radius_NDC * vec4(x / aspect, y, 0, 0)) * gl_in[0].gl_Position.w;
            gs_out.color = gs_in[0].color;                                     
            EmitVertex();                                               
        }

        void main() {    
            _position = gl_in[0].gl_Position / gl_in[0].gl_Position.w;
            emit(-1, -1);
            emit(1, -1);
            emit(-1, 1);
            emit(1, 1);
            EndPrimitive();
        }  
    )";

    char *_soup_frag_POINTS = R"(
        #version 330 core

        in GS_OUT {
            vec4 color;
            vec2 xy;
        } fs_in;

        out vec4 frag_color;

        void main() {
            frag_color.rgb = fs_in.color.rgb;
            frag_color.a = 1.0;
            if (length(fs_in.xy) > 1) { discard; }
        }
    )";

    char *_soup_geom_LINES = R"(
        #version 330 core
        layout (lines) in;
        layout (triangle_strip, max_vertices = 4) out;
        uniform float aspect;
        uniform float primitive_radius_NDC;

        in BLOCK {
            vec4 color;
        } gs_in[];

        out BLOCK {
            vec4 color;
        } gs_out;

        void main() {    
            vec4 s = gl_in[0].gl_Position / gl_in[0].gl_Position.w;
            vec4 t = gl_in[1].gl_Position / gl_in[1].gl_Position.w;
            vec4 color_s = gs_in[0].color;
            vec4 color_t = gs_in[1].color;

            vec4 perp = vec4(primitive_radius_NDC * vec2(1 / aspect, 1) * normalize(vec2(-1 / aspect, 1) * (t - s).yx), 0, 0);

            gl_Position = (s - perp) * gl_in[0].gl_Position.w;
            gs_out.color = color_s;
            EmitVertex();

            gl_Position = (t - perp) * gl_in[1].gl_Position.w;
            gs_out.color = color_t;
            EmitVertex();

            gl_Position = (s + perp) * gl_in[0].gl_Position.w;
            gs_out.color = color_s;
            EmitVertex();

            gl_Position = (t + perp) * gl_in[1].gl_Position.w;
            gs_out.color = color_t;
            EmitVertex();

            EndPrimitive();
        }  
    )";

    char *_soup_frag = R"(
        #version 330 core

        in BLOCK {
            vec4 color;
        } fs_in;

        out vec4 frag_color;

        void main() {
            frag_color = fs_in.color;
        }
    )";

    char *_itri_vert = R"(
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
        uniform vec3 color_if_vertex_colors_is_NULL;

        uniform bool has_vertex_texCoords;

        uniform mat4 P, V, M;

        void main() {
            vec4 tmp = M * vec4(vertex, 1);
            vs_out.position_World = vec3(tmp);
            gl_Position = P * V * tmp;
            vs_out.normal_World = inverse(transpose(mat3(M))) * normal;
            vs_out.color = has_vertex_colors ? color : color_if_vertex_colors_is_NULL;
            vs_out.texCoord = has_vertex_texCoords ? texCoord : vec2(0., 0.);
        }
    )";

    char *_itri_frag = R"(
        #version 330 core

        in BLOCK {
            vec3 position_World;
            vec3 normal_World;
            vec3 color;
            vec2 texCoord;
        } fs_in;
        uniform sampler2D i_texture;

        uniform vec4 eye_World;

        uniform bool has_vertex_normals;
        uniform bool has_vertex_colors;
        uniform bool has_vertex_texCoords;
        uniform bool has_texture;

        out vec4 frag_color;

        void main() {
            vec3 world_to_eye = vec3(eye_World) - fs_in.position_World;
            vec3 N = normalize(fs_in.normal_World);
            vec3 E = normalize(world_to_eye);

            vec3 rgb = fs_in.color;
            float a = 1.;

            if (has_texture) {

                vec2 texCoords = (has_vertex_texCoords) ? fs_in.texCoord : (.5 + .5 * N).xy;
                vec4 rgba = texture(i_texture, texCoords);
                rgb = rgba.rgb;
                a = rgba.a;

            } else if (has_vertex_normals) {

                rgb *= .8;

                float distance = length(world_to_eye);
                float attenuation = 1 / (1 + .02 * distance + .002 * distance * distance);

                vec3 L = E;
                vec3 H = normalize(L + E);
                float F0 = .05;

                float diffuse = max(0, dot(N, L));
                float specular = pow(max(0, dot(N, H)), 100);
                float fresnel = F0 + (1 - F0) * pow(1 - max(0, dot(N, H)), 5);

                rgb += attenuation * .7 * diffuse * vec3(.7, 1, .7);
                rgb += attenuation * .7 * specular * vec3(1, .2, .2);
                rgb += attenuation * .8 * fresnel * vec3(.2, .2, 1);

            }

            frag_color = vec4(rgb, a);
        }
    )";
};

struct C0_PersistsAcrossApps_NeverAutomaticallyClearedToZero__ManageItYourself {
    bool _app__completedFirstPass_helper;
    bool _app_completedFirstPass;
    bool _app_menu;
    int  _app_numApps;
    int  _app_index;
    int  _app_check;
    char _app_buffer[64];

    bool _cow_initialized;

    real *_eso_vertex_positions;
    real *_eso_vertex_colors;

    int _itri_shader_program;
    u32 _itri_VAO;
    u32 _itri_VBO[4];
    u32 _itri_EBO;

    FILE *_recorder_fp_mpg;
    FILE *_recorder_fp_dat;
    bool  _recorder_recording;
    int   _recorder_width;
    int   _recorder_height;
    int   _recorder_size_of_frame;
    u8   *_recorder__buffer1;
    u8   *_recorder_buffer2;
    int   _recorder_num_frames_recorded;
    bool  _recorder_out_of_space;

    int _soup_shader_program_POINTS;
    int _soup_shader_program_LINES;
    int _soup_shader_program_TRIANGLES;
    u32 _soup_VAO[3];
    u32 _soup_VBO[2];
    u32 _soup_EBO[3];

    GLFWwindow *_window_glfw_window;
    void       *_window_hwnd__note_this_is_NULL_if_not_on_Windows;
    real        _window_macbook_retina_scale;
};

struct C1_PersistsAcrossFrames_AutomaticallyClearedToZeroBetweenAppsBycow_reset {
    bool  _eso_called_eso_begin_before_calling_eso_vertex_or_eso_end;
    bool  _eso_called_eso_color_at_least_once_before_calling_eso_vertex;
    real  _eso_current_color[4];
    real  _eso_PVM[16];
    int   _eso_primitive;
    int   _eso_num_vertices;
    real  _eso_size_in_pixels;
    bool  _eso_use_world_units_instead_of_pixels;
    bool  _eso_overlay;

    real  _gui_x_curr;
    real  _gui_y_curr;
    void *_gui_selected;
    void *_gui_hot;
    real  _gui__dx_accumulator;
    bool  _gui_hide_and_disable;

    char _itri_texture_filenames[ITRI_MAX_NUM_TEXTURES][ITRI_MAX_FILENAME_LENGTH];
    u32  _itri_textures[ITRI_MAX_NUM_TEXTURES];
    int  _itri_num_textures;

    cs_audio_source_t *_sound_audio_source_ptrs[SOUND_MAX_DIFFERENT_FILES];
    char               _sound_filenames[SOUND_MAX_DIFFERENT_FILES][SOUND_MAX_FILENAME_LENGTH];
    int                _sound_num_loaded;
    real               _sound_music_gui_1_minus_volume;

    real _window_clear_color[4];
};


CW_USER_FACING_CONFIG config;
C2_READONLY_USER_FACING_DATA cow; // COW2
CX_INTERNAL_CONSTANTS COWX;
C0_PersistsAcrossApps_NeverAutomaticallyClearedToZero__ManageItYourself COW0;
C1_PersistsAcrossFrames_AutomaticallyClearedToZeroBetweenAppsBycow_reset COW1;


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// core ////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// #include "defines.cpp"///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// // assert
#define ASSERT(b) do { if (!(b)) { \
    printf("ASSERT Line %d in %s\n", __LINE__, __FILE__); \
    if (!config.tweaks_ASSERT_crashes_the_program_without_you_having_to_press_Enter) { \
        printf("press Enter to crash"); getchar(); \
    } \
    *((volatile int *) 0) = 0; \
} } while (0)
#define STATIC_ASSERT(cond) static_assert(cond, "STATIC_ASSERT");

// // working with real's
#define PI 3.14159265359
#define TINY_VAL 1e-7
#undef HUGE_VAL
#define HUGE_VAL 1e7
#define RAD(degrees) (PI / 180 * (degrees))
#define DEG(radians) (180. / PI * (radians))
#define INCHES(millimeters) ((millimeters) / 25.4)
#define MM(inches) ((inches) * 25.4)

#define ABS(a) ((a) < 0 ? -(a) : (a))
#define SGN(a) ((a) < 0 ? -1 : 1)
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define AVG(a, b) (.5 * (a) + .5 * (b))

#define IS_ZERO(a) (ABS(a) < TINY_VAL)
#define ARE_EQUAL(a, b) (IS_ZERO(ABS((a) - (b))))
#define IS_BETWEEN(p, a, b) (((a) - TINY_VAL < (p)) && ((p) < (b) + TINY_VAL))
#define IS_POSITIVE(a) ((a) > TINY_VAL)
#define IS_NEGATIVE(a) ((a) < -TINY_VAL)

#define LERP(t, a, b) ((1.0 - (t)) * (a) + (t) * (b)) // works on vecX, matX
#define INVERSE_LERP(p, a, b) (((p) - (a)) / real((b) - (a)))
#define LINEAR_REMAP(p, a, b, c, d) LERP(INVERSE_LERP(p, a, b), c, d)
#define CLAMP(t, a, b) MIN(MAX(t, a), b)
#define CLAMPED_LERP(t, a, b) LERP(CLAMP(t, 0.0, 1.0), a, b)
#define WRAP(t, a, b) ((a) + fmod((t) - (a), (b) - (a)))

// // working with int's
#define MODULO(x, N) (((x) % (N) + (N)) % (N)) // works on negative numbers
#define NELEMS(fixed_size_array) int(sizeof(fixed_size_array) / sizeof((fixed_size_array)[0]))
#define IS_ODD(a) ((a) % 2 != 0)
#define IS_EVEN(a) ((a) % 2 == 0)

// // misc.
#define STR(foo) #foo
#define XSTR(foo) STR(foo)
#define CONCAT_(a, b) a ## b
#define CONCAT(a, b) CONCAT_(a, b)
#define _NOOP(foo) foo

#define SUPPRESS_UNUSED_VARIABLE_COMPILER_WARNING(expr) do { (void)(expr); } while (0)

// // horrifying
#define do_once \
    static bool CONCAT(_do_once_, __LINE__) = false; \
bool CONCAT(_prev_do_once_, __LINE__) = CONCAT(_do_once_, __LINE__); \
CONCAT(_do_once_, __LINE__) = true; \
if (!CONCAT(_prev_do_once_, __LINE__) && CONCAT(_do_once_, __LINE__))

// // cmuratori filesizes
#define Kilobytes(Value) (          Value *1024LL)
#define Megabytes(Value) (Kilobytes(Value)*1024LL)
#define Gigabytes(Value) (Megabytes(Value)*1024LL)
#define Terabytes(Value) (Gigabytes(Value)*1024LL)

////////////////////////////////////////////////////////////////////////////////
// #include "_linalg.cpp"///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#define _LINALG_4X4(A, i, j) A[4 * (i) + (j)]

void _linalg_vec3_cross(real *c, real *a, real *b) { // c = a x b
    real tmp[3] = {
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    };
    memcpy(c, tmp, sizeof(tmp));
}

real _linalg_vecX_squared_length(int D, real *a) {
    real ret = 0;
    for (int d = 0; d < D; ++d) ret += pow(a[d], 2);
    return ret;
}

real _linalg_vecX_squared_distance(int D, real *a, real *b) {
    real ret = 0;
    for (int d = 0; d < D; ++d) ret += pow(a[d] - b[d], 2);
    return ret;
}

void _linalg_vecX_normalize(int D, real *a_hat, real *a) { // a_hat = a / |a|
    real L = sqrt(_linalg_vecX_squared_length(D, a));
    for (int d = 0; d < D; ++d) a_hat[d] = a[d] / L;
}

void _linalg_mat4_times_mat4(real *C, real *A, real *B) { // C = A B
    ASSERT(C);
    ASSERT(A);
    ASSERT(B);

    real tmp[16] = {}; { // allows for e.g. A <- A B;
        for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) for (int k = 0; k < 4; ++k) _LINALG_4X4(tmp, i, j) += _LINALG_4X4(A, i, k) * _LINALG_4X4(B, k, j);
    }
    memcpy(C, tmp, sizeof(tmp));
}

void _linalg_mat4_times_vec4_persp_divide(real *b, real *A, real *x) { // b = A x
    ASSERT(b);
    ASSERT(A);
    ASSERT(x);

    real tmp[4] = {}; { // allows for x <- A x
        for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) tmp[i] += _LINALG_4X4(A, i, j) * x[j];
        if (!IS_ZERO(tmp[3])) {
            ASSERT(!IS_ZERO(x[3]));
            for (int i = 0; i < 4; ++i) tmp[i] /= tmp[3]; 
        }
    }
    memcpy(b, tmp, sizeof(tmp));
}

real _linalg_mat4_determinant(real *A) {
    ASSERT(A);
    real A2323 = _LINALG_4X4(A, 2, 2) * _LINALG_4X4(A, 3, 3) - _LINALG_4X4(A, 2, 3) * _LINALG_4X4(A, 3, 2);
    real A1323 = _LINALG_4X4(A, 2, 1) * _LINALG_4X4(A, 3, 3) - _LINALG_4X4(A, 2, 3) * _LINALG_4X4(A, 3, 1);
    real A1223 = _LINALG_4X4(A, 2, 1) * _LINALG_4X4(A, 3, 2) - _LINALG_4X4(A, 2, 2) * _LINALG_4X4(A, 3, 1);
    real A0323 = _LINALG_4X4(A, 2, 0) * _LINALG_4X4(A, 3, 3) - _LINALG_4X4(A, 2, 3) * _LINALG_4X4(A, 3, 0);
    real A0223 = _LINALG_4X4(A, 2, 0) * _LINALG_4X4(A, 3, 2) - _LINALG_4X4(A, 2, 2) * _LINALG_4X4(A, 3, 0);
    real A0123 = _LINALG_4X4(A, 2, 0) * _LINALG_4X4(A, 3, 1) - _LINALG_4X4(A, 2, 1) * _LINALG_4X4(A, 3, 0);
    return _LINALG_4X4(A, 0, 0) * (_LINALG_4X4(A, 1, 1) * A2323 - _LINALG_4X4(A, 1, 2) * A1323 + _LINALG_4X4(A, 1, 3) * A1223) 
        -  _LINALG_4X4(A, 0, 1) * (_LINALG_4X4(A, 1, 0) * A2323 - _LINALG_4X4(A, 1, 2) * A0323 + _LINALG_4X4(A, 1, 3) * A0223) 
        +  _LINALG_4X4(A, 0, 2) * (_LINALG_4X4(A, 1, 0) * A1323 - _LINALG_4X4(A, 1, 1) * A0323 + _LINALG_4X4(A, 1, 3) * A0123) 
        -  _LINALG_4X4(A, 0, 3) * (_LINALG_4X4(A, 1, 0) * A1223 - _LINALG_4X4(A, 1, 1) * A0223 + _LINALG_4X4(A, 1, 2) * A0123);
}

void _linalg_mat4_inverse(real *invA, real *A) {
    ASSERT(invA);
    ASSERT(A);
    real one_over_det = 1 / _linalg_mat4_determinant(A);
    real tmp[16] = {}; { // allows for A <- inv(A)
        real A2323 = _LINALG_4X4(A, 2, 2) * _LINALG_4X4(A, 3, 3) - _LINALG_4X4(A, 2, 3) * _LINALG_4X4(A, 3, 2);
        real A1323 = _LINALG_4X4(A, 2, 1) * _LINALG_4X4(A, 3, 3) - _LINALG_4X4(A, 2, 3) * _LINALG_4X4(A, 3, 1);
        real A1223 = _LINALG_4X4(A, 2, 1) * _LINALG_4X4(A, 3, 2) - _LINALG_4X4(A, 2, 2) * _LINALG_4X4(A, 3, 1);
        real A0323 = _LINALG_4X4(A, 2, 0) * _LINALG_4X4(A, 3, 3) - _LINALG_4X4(A, 2, 3) * _LINALG_4X4(A, 3, 0);
        real A0223 = _LINALG_4X4(A, 2, 0) * _LINALG_4X4(A, 3, 2) - _LINALG_4X4(A, 2, 2) * _LINALG_4X4(A, 3, 0);
        real A0123 = _LINALG_4X4(A, 2, 0) * _LINALG_4X4(A, 3, 1) - _LINALG_4X4(A, 2, 1) * _LINALG_4X4(A, 3, 0);
        real A2313 = _LINALG_4X4(A, 1, 2) * _LINALG_4X4(A, 3, 3) - _LINALG_4X4(A, 1, 3) * _LINALG_4X4(A, 3, 2);
        real A1313 = _LINALG_4X4(A, 1, 1) * _LINALG_4X4(A, 3, 3) - _LINALG_4X4(A, 1, 3) * _LINALG_4X4(A, 3, 1);
        real A1213 = _LINALG_4X4(A, 1, 1) * _LINALG_4X4(A, 3, 2) - _LINALG_4X4(A, 1, 2) * _LINALG_4X4(A, 3, 1);
        real A2312 = _LINALG_4X4(A, 1, 2) * _LINALG_4X4(A, 2, 3) - _LINALG_4X4(A, 1, 3) * _LINALG_4X4(A, 2, 2);
        real A1312 = _LINALG_4X4(A, 1, 1) * _LINALG_4X4(A, 2, 3) - _LINALG_4X4(A, 1, 3) * _LINALG_4X4(A, 2, 1);
        real A1212 = _LINALG_4X4(A, 1, 1) * _LINALG_4X4(A, 2, 2) - _LINALG_4X4(A, 1, 2) * _LINALG_4X4(A, 2, 1);
        real A0313 = _LINALG_4X4(A, 1, 0) * _LINALG_4X4(A, 3, 3) - _LINALG_4X4(A, 1, 3) * _LINALG_4X4(A, 3, 0);
        real A0213 = _LINALG_4X4(A, 1, 0) * _LINALG_4X4(A, 3, 2) - _LINALG_4X4(A, 1, 2) * _LINALG_4X4(A, 3, 0);
        real A0312 = _LINALG_4X4(A, 1, 0) * _LINALG_4X4(A, 2, 3) - _LINALG_4X4(A, 1, 3) * _LINALG_4X4(A, 2, 0);
        real A0212 = _LINALG_4X4(A, 1, 0) * _LINALG_4X4(A, 2, 2) - _LINALG_4X4(A, 1, 2) * _LINALG_4X4(A, 2, 0);
        real A0113 = _LINALG_4X4(A, 1, 0) * _LINALG_4X4(A, 3, 1) - _LINALG_4X4(A, 1, 1) * _LINALG_4X4(A, 3, 0);
        real A0112 = _LINALG_4X4(A, 1, 0) * _LINALG_4X4(A, 2, 1) - _LINALG_4X4(A, 1, 1) * _LINALG_4X4(A, 2, 0);

        int i = 0;
        tmp[i++] =  one_over_det * (_LINALG_4X4(A, 1, 1) * A2323 - _LINALG_4X4(A, 1, 2) * A1323 + _LINALG_4X4(A, 1, 3) * A1223);
        tmp[i++] = -one_over_det * (_LINALG_4X4(A, 0, 1) * A2323 - _LINALG_4X4(A, 0, 2) * A1323 + _LINALG_4X4(A, 0, 3) * A1223);
        tmp[i++] =  one_over_det * (_LINALG_4X4(A, 0, 1) * A2313 - _LINALG_4X4(A, 0, 2) * A1313 + _LINALG_4X4(A, 0, 3) * A1213);
        tmp[i++] = -one_over_det * (_LINALG_4X4(A, 0, 1) * A2312 - _LINALG_4X4(A, 0, 2) * A1312 + _LINALG_4X4(A, 0, 3) * A1212);
        tmp[i++] = -one_over_det * (_LINALG_4X4(A, 1, 0) * A2323 - _LINALG_4X4(A, 1, 2) * A0323 + _LINALG_4X4(A, 1, 3) * A0223);
        tmp[i++] =  one_over_det * (_LINALG_4X4(A, 0, 0) * A2323 - _LINALG_4X4(A, 0, 2) * A0323 + _LINALG_4X4(A, 0, 3) * A0223);
        tmp[i++] = -one_over_det * (_LINALG_4X4(A, 0, 0) * A2313 - _LINALG_4X4(A, 0, 2) * A0313 + _LINALG_4X4(A, 0, 3) * A0213);
        tmp[i++] =  one_over_det * (_LINALG_4X4(A, 0, 0) * A2312 - _LINALG_4X4(A, 0, 2) * A0312 + _LINALG_4X4(A, 0, 3) * A0212);
        tmp[i++] =  one_over_det * (_LINALG_4X4(A, 1, 0) * A1323 - _LINALG_4X4(A, 1, 1) * A0323 + _LINALG_4X4(A, 1, 3) * A0123);
        tmp[i++] = -one_over_det * (_LINALG_4X4(A, 0, 0) * A1323 - _LINALG_4X4(A, 0, 1) * A0323 + _LINALG_4X4(A, 0, 3) * A0123);
        tmp[i++] =  one_over_det * (_LINALG_4X4(A, 0, 0) * A1313 - _LINALG_4X4(A, 0, 1) * A0313 + _LINALG_4X4(A, 0, 3) * A0113);
        tmp[i++] = -one_over_det * (_LINALG_4X4(A, 0, 0) * A1312 - _LINALG_4X4(A, 0, 1) * A0312 + _LINALG_4X4(A, 0, 3) * A0112);
        tmp[i++] = -one_over_det * (_LINALG_4X4(A, 1, 0) * A1223 - _LINALG_4X4(A, 1, 1) * A0223 + _LINALG_4X4(A, 1, 2) * A0123);
        tmp[i++] =  one_over_det * (_LINALG_4X4(A, 0, 0) * A1223 - _LINALG_4X4(A, 0, 1) * A0223 + _LINALG_4X4(A, 0, 2) * A0123);
        tmp[i++] = -one_over_det * (_LINALG_4X4(A, 0, 0) * A1213 - _LINALG_4X4(A, 0, 1) * A0213 + _LINALG_4X4(A, 0, 2) * A0113);
        tmp[i++] =  one_over_det * (_LINALG_4X4(A, 0, 0) * A1212 - _LINALG_4X4(A, 0, 1) * A0212 + _LINALG_4X4(A, 0, 2) * A0112);
    }
    memcpy(invA, tmp, sizeof(tmp));
}

void _linalg_mat4_transpose(real *AT, real *A) { // AT = A^T
    ASSERT(AT);
    ASSERT(A);
    real tmp[16] = {}; { // allows for A <- transpose(A)
        for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) _LINALG_4X4(tmp, i, j) = _LINALG_4X4(A, j, i);
    }
    memcpy(AT, tmp, 16 * sizeof(real));
}

////////////////////////////////////////////////////////////////////////////////
// #include "window.cpp"////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void _window_swap_draw_buffers() {
    ASSERT(COW0._window_glfw_window);
    glfwSwapBuffers(COW0._window_glfw_window);
}

void _window_clear_draw_buffer() {
    ASSERT(COW0._window_glfw_window);
    glClearColor(
            float(COW1._window_clear_color[0]),
            float(COW1._window_clear_color[1]),
            float(COW1._window_clear_color[2]),
            float(COW1._window_clear_color[3])
            );
    glClear(
            GL_COLOR_BUFFER_BIT
            | GL_DEPTH_BUFFER_BIT
            | GL_STENCIL_BUFFER_BIT
           );
}

void window_set_title(char *title) {
    ASSERT(COW0._window_glfw_window);
    ASSERT(title);
    glfwSetWindowTitle(COW0._window_glfw_window, title);
}

void window_set_position(real x, real y) {
    ASSERT(COW0._window_glfw_window);
    glfwSetWindowPos(COW0._window_glfw_window, int(x), int(y));
}

void window_set_size(real width, real height) {
    ASSERT(COW0._window_glfw_window);
    glfwSetWindowSize(COW0._window_glfw_window, int(width), int(height));
}

void window_set_height__16_by_9_aspect(real height) {
    window_set_size(height * 16 / 9, height);
}

void window_set_clear_color(real r, real g, real b, real a = 1.0) {
    COW1._window_clear_color[0] = r;
    COW1._window_clear_color[1] = g;
    COW1._window_clear_color[2] = b;
    COW1._window_clear_color[3] = a;
}

void window_set_floating(bool floating) {
    glfwSetWindowAttrib(COW0._window_glfw_window, GLFW_FLOATING, floating);
}

void window_set_decorated(bool decorated) {
    glfwSetWindowAttrib(COW0._window_glfw_window, GLFW_DECORATED, decorated);
}

void _callback_set_callbacks();
void _window_init() {
    ASSERT(glfwInit());

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    // glfwWindowHint(GLFW_TRANSPARENT_FRAMEBUFFER, GLFW_TRUE);

    // glfwWindowHint(GLFW_SAMPLES, 4); // multisampling

    glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);

    COW0._window_glfw_window = glfwCreateWindow(1000, 1000, "", NULL, NULL);
    if (!COW0._window_glfw_window) {
        printf("[cow] something's gone wonky; if you weren't just messing with init(...) or something, please try restarting your computer and try again.\n");
        ASSERT(0);
    }

    {
        GLFWimage images[1]; 
        images[0].pixels = stbi_load("codebase/icon.png", &images[0].width, &images[0].height, 0, 4); //rgba channels 
        if (images[0].pixels) {
            glfwSetWindowIcon(COW0._window_glfw_window, 1, images); 
            stbi_image_free(images[0].pixels);
        }
    }

    glfwMakeContextCurrent(COW0._window_glfw_window);

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

    glDisable(GL_CULL_FACE);   // fornow

    glfwSwapInterval(1);

    _callback_set_callbacks();

    #if defined(WIN32) || defined(_WIN64)
    COW0._window_hwnd__note_this_is_NULL_if_not_on_Windows = glfwGetWin32Window(COW0._window_glfw_window);
    #endif

    { // _macbook_retina_scale
        int num, den, _;
        glfwGetFramebufferSize(COW0._window_glfw_window, &num, &_);
        glfwGetWindowSize(COW0._window_glfw_window, &den, &_);
        COW0._window_macbook_retina_scale = num / den;
    }
}

void _window_reset() {
    glfwSetInputMode(COW0._window_glfw_window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
    COW1._window_clear_color[0] = config._default_window_color_rgb[0];
    COW1._window_clear_color[1] = config._default_window_color_rgb[1];
    COW1._window_clear_color[2] = config._default_window_color_rgb[2];
    COW1._window_clear_color[3] = config.default_window_color_a;
    window_set_position(
            config.default_window_position[0],
            config.default_window_position[1]
            );
    window_set_size(
            config.default_window_size[0],
            config.default_window_size[1]
            );
    window_set_floating(config.default_window_floating);
    window_set_decorated(config.default_window_decorated);
    _window_clear_draw_buffer();
    glfwShowWindow(COW0._window_glfw_window); // only actually does anything on first reset
}

void _window_begin_frame() {
    _window_swap_draw_buffers();
    _window_clear_draw_buffer();
}

void _window_get_size(real *width, real *height) {
    ASSERT(COW0._window_glfw_window);
    int _width, _height;
    glfwGetFramebufferSize(COW0._window_glfw_window, &_width, &_height);
    *width = real(_width);
    *height = real(_height);
}
real window_get_height_in_pixels() {
    ASSERT(COW0._window_glfw_window);
    real _, height;
    _window_get_size(&_, &height);
    return height;
}
real window_get_aspect() {
    real width, height;
    _window_get_size(&width, &height);
    return width / height;
}

#ifdef SNAIL_CPP
vec2 window_get_size() {
    real width, height;
    _window_get_size(&width, &height);
    return { width, height };
}

void window_set_size(vec2 size) {
    window_set_size(size[0], size[1]);
}

void window_set_position(vec2 position) {
    window_set_position(position[0], position[1]);
}

void window_set_clear_color(vec3 rgb, real a = 1.0) {
    COW1._window_clear_color[0] = rgb[0];
    COW1._window_clear_color[1] = rgb[1];
    COW1._window_clear_color[2] = rgb[2];
    COW1._window_clear_color[3] = a;
}

#endif

////////////////////////////////////////////////////////////////////////////////
// #include "window_tform.cpp"//////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void _window_get_P_perspective(real *P, real angle_of_view, real n = 0, real f = 0, real aspect = 0) {
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

    real angle_y = angle_of_view / 2;
    real Q_y = 1 / tan(angle_y);
    real Q_x = Q_y / aspect;

    memset(P, 0, 16 * sizeof(real));
    _LINALG_4X4(P, 0, 0) = Q_x;
    _LINALG_4X4(P, 0, 1) = 0; // TERRIBLE PATCH
    _LINALG_4X4(P, 1, 1) = Q_y;
    _LINALG_4X4(P, 3, 2) = -1;

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
    _LINALG_4X4(P, 2, 2) = (n + f) / (n - f);     // a
    _LINALG_4X4(P, 2, 3) = (2 * n * f) / (f - n); // b
}

void _window_get_P_ortho(real *P, real screen_height_World, real n = 0, real f = 0, real aspect = 0) {
    ASSERT(P);
    ASSERT(!IS_ZERO(screen_height_World));
    if (ARE_EQUAL(n, f)) {
        n = 1000.0;
        f =  -1000.0; 
    }
    if (IS_ZERO(aspect)) { aspect = window_get_aspect(); }

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

    real r_y = screen_height_World / 2;
    real r_x = window_get_aspect() * r_y;

    // [x'] = [1/r_x      0   0  0] [x] = [ x/r_x]
    // [y'] = [    0  1/r_y   0  0] [y] = [ y/r_y]
    // [z'] = [    0      0   a  b] [z] = [az + b]
    // [1 ] = [    0      0   0  1] [1] = [     1]

    memset(P, 0, 16 * sizeof(real));
    _LINALG_4X4(P, 0, 0) = 1 / r_x;
    _LINALG_4X4(P, 1, 1) = 1 / r_y;

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
    _LINALG_4X4(P, 2, 2) = 2 / (f - n);
    _LINALG_4X4(P, 2, 3) = (n + f) / (n - f);

    // ???
    // _LINALG_4X4(P, 2, 2) *= -1; // sign flip
    // _LINALG_4X4(P, 2, 3) *= -1; //          
    _LINALG_4X4(P, 3, 3) = 1;
}

void _window_get_NDC_from_Screen(real *NDC_from_Screen) {
    _window_get_P_ortho(NDC_from_Screen, window_get_height_in_pixels());
    _LINALG_4X4(NDC_from_Screen, 1, 1) *= -1;
    _LINALG_4X4(NDC_from_Screen, 0, 3) -= 1;
    _LINALG_4X4(NDC_from_Screen, 1, 3) += 1;

    _LINALG_4X4(NDC_from_Screen, 0, 0) *= COW0._window_macbook_retina_scale;
    _LINALG_4X4(NDC_from_Screen, 1, 1) *= COW0._window_macbook_retina_scale;
}

void _window_get_PV_ortho(real *PV, real screen_height_World, real eye_x, real eye_y) {
    ASSERT(PV);
    // fornow slow
    real P[16];
    _window_get_P_ortho(P, screen_height_World);
    real V[16] = { 
        1, 0, 0, -eye_x,
        0, 1, 0, -eye_y,
        0, 0, 1, 0     ,
        0, 0, 0, 1     ,
    };
    _linalg_mat4_times_mat4(PV, P, V);
}

#ifdef SNAIL_CPP
mat4 _window_get_P_perspective(real angle_of_view, real n = 0, real f = 0, real aspect = 0) {
    mat4 ret;
    _window_get_P_perspective(ret.data, angle_of_view, n, f, aspect);
    return ret;
}

mat4 _window_get_P_ortho(real screen_height_World, real n = 0, real f = 0, real aspect = 0) {
    mat4 ret;
    _window_get_P_ortho(ret.data, screen_height_World, n, f, aspect);
    return ret;
}

mat4 _window_get_NDC_from_Screen() {
    mat4 ret;
    _window_get_NDC_from_Screen(ret.data);
    return ret;
}
#endif

////////////////////////////////////////////////////////////////////////////////
// #include "input_and_callback.cpp"////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void _input_begin_frame() {
    ASSERT(COW0._cow_initialized);
    memset(cow.key_pressed, 0, sizeof(cow.key_pressed));
    memset(cow.key_released, 0, sizeof(cow.key_released));
    cow.mouse_left_pressed = false;
    cow.mouse_left_released = false;
    cow.mouse_right_pressed = false;
    cow.mouse_right_released = false;
    cow.mouse_change_in_position_Screen[0] = 0;
    cow.mouse_change_in_position_Screen[1] = 0;
    cow.mouse_change_in_position_NDC[0] = 0;
    cow.mouse_change_in_position_NDC[1] = 0;
    cow.mouse_wheel_offset = 0;
    glfwPollEvents();
    // NOTE e.g. key_*['j'] the same as key_*['J']
    for (int i = 0; i < 26; ++i) {
        cow.key_pressed ['a' + i] = cow.key_pressed ['A' + i];
        cow.key_held    ['a' + i] = cow.key_held    ['A' + i];
        cow.key_released['a' + i] = cow.key_released['A' + i];
    }
}

void _input_get_mouse_position_and_change_in_position_in_world_coordinates(
        real *PV,
        real *mouse_x_world,
        real *mouse_y_world,
        real *mouse_dx_world,
        real *mouse_dy_world) {
    ASSERT(PV);
    real World_from_NDC[16] = {};
    _linalg_mat4_inverse(World_from_NDC, PV);
    real   xy[4] = { cow.mouse_position_NDC[0],  cow.mouse_position_NDC[1],  0.0, 1.0 }; // point
    real dxdy[4] = { cow.mouse_change_in_position_NDC[0], cow.mouse_change_in_position_NDC[1], 0.0, 0.0 }; // vector
    _linalg_mat4_times_vec4_persp_divide(xy, World_from_NDC, xy);
    _linalg_mat4_times_vec4_persp_divide(dxdy, World_from_NDC, dxdy);
    if (mouse_x_world) *mouse_x_world = xy[0];
    if (mouse_y_world) *mouse_y_world = xy[1];
    if (mouse_dx_world) *mouse_dx_world = dxdy[0];
    if (mouse_dy_world) *mouse_dy_world = dxdy[1];
}

#ifdef SNAIL_CPP
vec2 input_get_mouse_position_in_world_coordinates(mat4 PV) {
    vec2 ret;
    _input_get_mouse_position_and_change_in_position_in_world_coordinates(PV.data, &ret[0], &ret[1], NULL, NULL);
    return ret;
}
vec2 input_get_mouse_change_in_position_in_world_coordinates(mat4 PV) {
    vec2 ret;
    _input_get_mouse_position_and_change_in_position_in_world_coordinates(PV.data, NULL, NULL, &ret[0], &ret[1]);
    return ret;
}
#endif

void _callback_key(GLFWwindow *, int key, int, int action, int mods) {
    if (key < 0) { return; }
    if (key >= 512) { return; }
    if (action == GLFW_PRESS) {
        cow.key_pressed[key] = true;
        cow.key_held[key] = true;
        // cow.key_toggle[key] = !cow.key_toggle[key];
    } else if (action == GLFW_RELEASE) {
        cow.key_released[key] = true;
        cow.key_held[key] = false;
    }
    cow.key_shift_held = (mods & GLFW_MOD_SHIFT);
}

void _callback_cursor_position(GLFWwindow *, real xpos, real ypos) {
    real tmp_mouse_s_NDC_0 = cow.mouse_position_NDC[0];
    real tmp_mouse_s_NDC_1 = cow.mouse_position_NDC[1];
    real tmp_mouse_s_Screen_0 = cow.mouse_position_Screen[0];
    real tmp_mouse_s_Screen_1 = cow.mouse_position_Screen[1];

    { // macbook retina nonsense
        xpos *= COW0._window_macbook_retina_scale;
        ypos *= COW0._window_macbook_retina_scale;
    }

    cow.mouse_position_Screen[0] = xpos;
    cow.mouse_position_Screen[1] = ypos;
    real s_NDC[4] = {}; {
        real s_Screen[4] = { xpos, ypos, 0, 1 };
        _linalg_mat4_times_vec4_persp_divide(s_NDC, (real *) &cow.NDC_from_Screen, s_Screen);
    }
    cow.mouse_position_NDC[0] = s_NDC[0];
    cow.mouse_position_NDC[1] = s_NDC[1];
    // callback may be called multiple times per frame
    cow.mouse_change_in_position_Screen[0] += (cow.mouse_position_Screen[0] - tmp_mouse_s_Screen_0);
    cow.mouse_change_in_position_Screen[1] += (cow.mouse_position_Screen[1] - tmp_mouse_s_Screen_1);
    cow.mouse_change_in_position_NDC[0] += (cow.mouse_position_NDC[0] - tmp_mouse_s_NDC_0);
    cow.mouse_change_in_position_NDC[1] += (cow.mouse_position_NDC[1] - tmp_mouse_s_NDC_1);
}

void _callback_mouse_button(GLFWwindow *, int button, int action, int) {
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        if (action == GLFW_PRESS) { 
            cow.mouse_left_pressed = true;
            cow.mouse_left_held = true;
        } else if (action == GLFW_RELEASE) { 
            cow.mouse_left_released = true;
            cow.mouse_left_held = false;
        }
    } else if (button == GLFW_MOUSE_BUTTON_RIGHT) {
        if (action == GLFW_PRESS) { 
            cow.mouse_right_pressed = true;
            cow.mouse_right_held = true;
        } else if (action == GLFW_RELEASE) { 
            cow.mouse_right_held = false;
        }
    }
}

void _callback_scroll(GLFWwindow *, real, real yoffset) {
    cow.mouse_wheel_offset += yoffset;
}

void _callback_framebuffer_size(GLFWwindow *, int width, int height) {
    glViewport(0, 0, width, height);
}

void _callback_set_callbacks() {
    glfwSetFramebufferSizeCallback(COW0._window_glfw_window, _callback_framebuffer_size);
    glfwSetKeyCallback(COW0._window_glfw_window, _callback_key);
    glfwSetCursorPosCallback(COW0._window_glfw_window, _callback_cursor_position);
    glfwSetMouseButtonCallback(COW0._window_glfw_window, _callback_mouse_button);
    glfwSetScrollCallback(COW0._window_glfw_window, _callback_scroll);
}

////////////////////////////////////////////////////////////////////////////////
// #include "app.cpp"///////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

bool cow_begin_frame();
bool _app_while_loop_condition() {
    if (!COW0._app__completedFirstPass_helper) {
        COW0._app__completedFirstPass_helper = true;
    } else {
        COW0._app_completedFirstPass = true;
    }
    COW0._app_check = 0;
    cow_begin_frame();
    return true;
}

#define APPS \
    if (!COW0._cow_initialized) { _cow_init(); } \
    COW0._app__completedFirstPass_helper = false; \
    COW0._app_completedFirstPass = false; \
    COW0._app_menu = false; \
    COW0._app_numApps = 0; \
    COW0._app_index = 0; \
    while (_app_while_loop_condition())

#define APP(_app_name) \
    if (!COW0._app_completedFirstPass) { \
        ++COW0._app_numApps; \
    } else { \
        if (COW0._app_menu) { \
            if (cow.key_pressed['q']) { break; } \
            if (gui_button(STR(_app_name)"()", 'a' + COW0._app_check)) { \
                COW0._app_index = COW0._app_check; \
                COW0._app_menu = false; \
            } \
        } else { \
            if (COW0._app_index == COW0._app_check) { \
                sprintf(COW0._app_buffer, "%d/%d - %s", COW0._app_index + 1, COW0._app_numApps, STR(_app_name)); \
                _cow_reset(); \
                window_set_title(COW0._app_buffer); \
                _app_name(); \
                if (cow.key_pressed['q'] && cow.key_shift_held) {                                                          break;                                 } \
                if (cow.key_pressed['q'                 ]) { COW0._app_index++; if (COW0._app_index == COW0._app_numApps) { break;                               } } \
                if (cow.key_pressed[config.hotkeys_app_next]) { COW0._app_index++; if (COW0._app_index == COW0._app_numApps) { COW0._app_index = 0;                   } } \
                if (cow.key_pressed[config.hotkeys_app_prev]) { COW0._app_index--; if (COW0._app_index == -1             ) { COW0._app_index = COW0._app_numApps - 1; } } \
                if (cow.key_pressed[config.hotkeys_app_menu]) { \
                    COW0._app_index = 0; \
                    COW0._app_menu = true; \
                    _cow_reset(); \
                    window_set_title("menu"); \
                    continue; \
                } \
            } \
        } \
        ++COW0._app_check; \
    }

////////////////////////////////////////////////////////////////////////////////
// #include "shader.cpp"////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int _shader_compile(char *source, GLenum type) {
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

int _shader_build_program(int vertex_shader, int fragment_shader, int geometry_shader = 0) {
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

int shader_compile_and_build_program(char *vertex_shader_source, char *fragment_shader_source, char *geometry_shader_source = 0) {
    int vert = _shader_compile(vertex_shader_source, GL_VERTEX_SHADER);
    int frag = _shader_compile(fragment_shader_source, GL_FRAGMENT_SHADER);
    int geom = geometry_shader_source ? _shader_compile(geometry_shader_source, GL_GEOMETRY_SHADER) : 0;
    return _shader_build_program(vert, frag, geom);
}

int _shader_get_uniform_location(int ID, char *name) {
    ASSERT(ID);
    int location = glGetUniformLocation(ID, name);
    return location;
}

void _shader_set_uniform_real(int ID, char *name, real value) {
    glUniform1f(_shader_get_uniform_location(ID, name), float(value));
}

void _shader_set_uniform_int(int ID, char *name, int value) {
    glUniform1i(_shader_get_uniform_location(ID, name), value);
}

void _shader_set_uniform_bool(int ID, char *name, bool value) {
    glUniform1ui(_shader_get_uniform_location(ID, name), value);
}

void _shader_set_uniform_vec2(int ID, char *name, real *value) {
    ASSERT(value);
    glUniform2f(_shader_get_uniform_location(ID, name), float(value[0]), float(value[1]));
}

void _shader_set_uniform_vec3(int ID, char *name, real *value) {
    ASSERT(value);
    glUniform3f(_shader_get_uniform_location(ID, name), float(value[0]), float(value[1]), float(value[2]));
}

void _shader_set_uniform_vec4(int ID, char *name, real *value) {
    ASSERT(value);
    glUniform4f(_shader_get_uniform_location(ID, name), float(value[0]), float(value[1]), float(value[2]), float(value[3]));
}

void _shader_set_uniform_mat4(int ID, char *name, real *value) {
    ASSERT(value);
    float as_floats[16] = {}; {
        for (int k = 0; k < 16; ++k) as_floats[k] = float(value[k]);
    }
    glUniformMatrix4fv(_shader_get_uniform_location(ID, name), 1, GL_TRUE, as_floats);
}

void _shader_set_uniform_array_vec3(int ID, char *name, int count, real *value) {
    ASSERT(value);
    float *tmp = (float *) malloc(count * 3 * sizeof(float)); 
    for (int i = 0; i < count; ++i) {
        for (int d = 0; d < 3; ++d) { 
            tmp[3 * i + d] = float(value[3 * i + d]);
        }
    }
    glUniform3fv(_shader_get_uniform_location(ID, name), count, tmp);
    free(tmp);
}

#ifdef SNAIL_CPP
void shader_set_uniform(int ID, char *name, real value) {
    _shader_set_uniform_real(ID, name, value);
}

void shader_set_uniform(int ID, char *name, int value)  {
    _shader_set_uniform_int(ID, name, value);
}

void shader_set_uniform(int ID, char *name, bool value) {
    _shader_set_uniform_bool(ID, name, value);
}

void shader_set_uniform(int ID, char *name, vec2 value) {
    _shader_set_uniform_vec2(ID, name, (real *) &value);
}

void shader_set_uniform(int ID, char *name, vec3 value) {
    _shader_set_uniform_vec3(ID, name, (real *) &value);
}

void shader_set_uniform(int ID, char *name, vec4 value) {
    _shader_set_uniform_vec4(ID, name, (real *) &value);
}

void shader_set_uniform(int ID, char *name, mat4 value) {
    _shader_set_uniform_mat4(ID, name, value.data);
}

void shader_set_uniform(int ID, char *name, int count, vec3 *value) {
    _shader_set_uniform_array_vec3(ID, name, count, (real *) value);
}
#endif

////////////////////////////////////////////////////////////////////////////////
// #include "soup.cpp"//////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#define SOUP_POINTS         GL_POINTS
#define SOUP_LINES          GL_LINES
#define SOUP_LINE_STRIP     GL_LINE_STRIP
#define SOUP_LINE_LOOP      GL_LINE_LOOP
#define SOUP_TRIANGLES      GL_TRIANGLES
#define SOUP_TRIANGLE_FAN   GL_TRIANGLE_FAN
#define SOUP_TRIANGLE_STRIP GL_TRIANGLE_STRIP
#define SOUP_TRIANGLE_MESH  254
#define SOUP_QUADS          255
#define SOUP_QUAD_MESH      253


#define _SOUP_XY 2
#define _SOUP_XYZ 3
#define _SOUP_XYZW 4
#define _SOUP_RGB 3
#define _SOUP_RGBA 4

void _soup_init() {
    COW0._soup_shader_program_POINTS = shader_compile_and_build_program(COWX._soup_vert, COWX._soup_frag_POINTS, COWX._soup_geom_POINTS);
    COW0._soup_shader_program_LINES = shader_compile_and_build_program(COWX._soup_vert, COWX._soup_frag, COWX._soup_geom_LINES);
    COW0._soup_shader_program_TRIANGLES = shader_compile_and_build_program(COWX._soup_vert, COWX._soup_frag);
    glGenVertexArrays(NELEMS(COW0._soup_VAO), COW0._soup_VAO);
    glGenBuffers(NELEMS(COW0._soup_VBO), COW0._soup_VBO);
    glGenBuffers(NELEMS(COW0._soup_EBO), COW0._soup_EBO);
}

void _soup_draw(
        real *PVM,
        int primitive,
        int dimension_of_positions,
        int dimension_of_colors,
        int num_vertices,
        real *vertex_positions,
        real *vertex_colors,
        real r_if_vertex_colors_is_NULL,
        real g_if_vertex_colors_is_NULL,
        real b_if_vertex_colors_is_NULL,
        real a_if_vertex_colors_is_NULL,
        real size_in_pixels,
        bool use_world_units_instead_of_pixels,
        bool force_draw_on_top) {

    if (num_vertices == 0) { return; } // NOTE: num_vertices zero is valid input

    ASSERT(PVM);
    ASSERT(dimension_of_positions >= 1 && dimension_of_positions <= 4);
    if (vertex_colors) ASSERT(dimension_of_colors == 3 || dimension_of_colors == 4);
    ASSERT(dimension_of_colors >= 0);
    ASSERT(vertex_positions);

    { // scrub input
        if (IS_ZERO(size_in_pixels)) { size_in_pixels = config.tweaks_size_in_pixels_soup_draw_defaults_to_if_you_pass_0_for_size_in_pixels; }
        size_in_pixels *= COW0._window_macbook_retina_scale;
    }

    int mesh_special_case = 0; {
        if (primitive == SOUP_TRIANGLE_MESH || primitive == SOUP_QUAD_MESH) {
            _soup_draw(
                    PVM, primitive == SOUP_TRIANGLE_MESH ? SOUP_TRIANGLES : SOUP_QUADS,
                    dimension_of_positions,
                    dimension_of_colors,
                    num_vertices,
                    vertex_positions,
                    vertex_colors,
                    r_if_vertex_colors_is_NULL,
                    g_if_vertex_colors_is_NULL,
                    b_if_vertex_colors_is_NULL,
                    a_if_vertex_colors_is_NULL,
                    size_in_pixels,
                    use_world_units_instead_of_pixels,
                    force_draw_on_top);

            if (primitive == SOUP_TRIANGLE_MESH) {
                mesh_special_case = 1;
            } else {
                mesh_special_case = 2;
            }

            primitive = SOUP_LINES;
            vertex_colors = NULL;
            r_if_vertex_colors_is_NULL = 1.0;
            g_if_vertex_colors_is_NULL = 1.0;
            b_if_vertex_colors_is_NULL = 1.0;
            a_if_vertex_colors_is_NULL = 1.0;
        }
    }

    if (config.tweaks_soup_draw_with_rounded_corners_for_all_line_primitives && (primitive == SOUP_LINES || primitive == SOUP_LINE_STRIP || primitive == SOUP_LINE_LOOP)) {
        _soup_draw(
                PVM,
                SOUP_POINTS,
                dimension_of_positions,
                dimension_of_colors,
                num_vertices,
                vertex_positions,
                vertex_colors,
                r_if_vertex_colors_is_NULL,
                g_if_vertex_colors_is_NULL,
                b_if_vertex_colors_is_NULL,
                a_if_vertex_colors_is_NULL,
                size_in_pixels,
                use_world_units_instead_of_pixels,
                force_draw_on_top);
    }

    real color_if_vertex_colors_is_NULL[4] = { r_if_vertex_colors_is_NULL, g_if_vertex_colors_is_NULL, b_if_vertex_colors_is_NULL, a_if_vertex_colors_is_NULL };

    glBindVertexArray(COW0._soup_VAO[mesh_special_case]);
    int i_attrib = 0;
    auto guarded_push = [&](int buffer_size, void *array, int dim) {
        glDisableVertexAttribArray(i_attrib); // fornow
        if (array) {
            glBindBuffer(GL_ARRAY_BUFFER, COW0._soup_VBO[i_attrib]);
            glBufferData(GL_ARRAY_BUFFER, buffer_size, array, GL_DYNAMIC_DRAW);
            glVertexAttribPointer(i_attrib, dim, GL_REAL, 0, 0, NULL);
            glEnableVertexAttribArray(i_attrib);
        }
        ++i_attrib;
    };
    int vvv_size = int(num_vertices * dimension_of_positions * sizeof(real));
    int ccc_size = int(num_vertices * dimension_of_colors * sizeof(real));
    guarded_push(vvv_size, vertex_positions, dimension_of_positions);
    guarded_push(ccc_size, vertex_colors, dimension_of_colors);

    int shader_program = 0; {
        if (primitive == SOUP_POINTS) {
            shader_program = COW0._soup_shader_program_POINTS;
        } else if (primitive == SOUP_LINES || primitive == SOUP_LINE_STRIP || primitive == SOUP_LINE_LOOP) {
            shader_program = COW0._soup_shader_program_LINES;
        } else { // including SOUP_QUADS
            shader_program = COW0._soup_shader_program_TRIANGLES;
        }
    }
    ASSERT(shader_program);
    glUseProgram(shader_program);

    _shader_set_uniform_real(shader_program, "aspect", window_get_aspect());
    if (use_world_units_instead_of_pixels) {
        real tmp[4] = {  0.0, 0.5 * size_in_pixels, 0.0 , 0.0 };
        _linalg_mat4_times_vec4_persp_divide(tmp, PVM, tmp);
        _shader_set_uniform_real(shader_program, "primitive_radius_NDC",
                sqrt(_linalg_vecX_squared_length(4, tmp)));
    } else {
        _shader_set_uniform_real(shader_program, "primitive_radius_NDC",
                0.5 * size_in_pixels / window_get_height_in_pixels());
    } 
    _shader_set_uniform_bool(shader_program, "has_vertex_colors", vertex_colors != NULL);
    _shader_set_uniform_bool(shader_program, "force_draw_on_top", force_draw_on_top);
    _shader_set_uniform_mat4(shader_program, "PVM", PVM);
    _shader_set_uniform_vec4(shader_program, "color_if_vertex_colors_is_NULL", color_if_vertex_colors_is_NULL);

    // glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    if (primitive != SOUP_QUADS && !mesh_special_case) {
        glDrawArrays(primitive, 0, num_vertices);
    } else {
        // we upload three EBO's _once_               
        // and bind the appropriate one before drawing

        ASSERT(primitive == SOUP_QUADS || mesh_special_case != 0);

        const int MAX_VERTICES = 1000000;
        ASSERT(num_vertices <= MAX_VERTICES);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, COW0._soup_EBO[mesh_special_case]);

        if (primitive == SOUP_QUADS) {
            primitive = SOUP_TRIANGLES;
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

#ifdef SNAIL_CPP
template <int D_pos, int D_color = 3> void soup_draw(
        mat4 PVM,
        int primitive,
        int num_vertices,
        Vec<D_pos> *vertex_positions,
        Vec<D_color> *vertex_colors,
        Vec<D_color> color_if_vertex_colors_is_NULL,
        real size_in_pixels = 0,
        bool use_world_units_instead_of_pixels = false,
        bool force_draw_on_top = false) {
    STATIC_ASSERT(D_pos == 2 || D_pos == 3 || D_pos == 4);
    STATIC_ASSERT(D_color == 3 || D_color == 4);

    _soup_draw(
            PVM.data,
            primitive,
            D_pos,
            D_color,
            num_vertices,
            (real *) vertex_positions,
            (real *) vertex_colors,
            color_if_vertex_colors_is_NULL[0],
            color_if_vertex_colors_is_NULL[1],
            color_if_vertex_colors_is_NULL[2],
            (D_color == 4) ? color_if_vertex_colors_is_NULL[3] : 1,
            size_in_pixels,
            use_world_units_instead_of_pixels,
            force_draw_on_top
            );
}

template <int D_pos, int D_color> void soup_draw(
        mat4 PVM,
        int primitive,
        int num_vertices,
        Vec<D_pos> *vertex_positions,
        void *vertex_colors,
        Vec<D_color> color_if_vertex_colors_is_NULL,
        real size_in_pixels = 0,
        bool use_world_units_instead_of_pixels = false,
        bool force_draw_on_top = false) {

    ASSERT(vertex_colors == NULL);

    _soup_draw(
            PVM.data,
            primitive,
            D_pos,
            D_color,
            num_vertices,
            (real *) vertex_positions,
            (real *) vertex_colors,
            color_if_vertex_colors_is_NULL[0],
            color_if_vertex_colors_is_NULL[1],
            color_if_vertex_colors_is_NULL[2],
            (D_color == 4) ? color_if_vertex_colors_is_NULL[3] : 1,
            size_in_pixels,
            use_world_units_instead_of_pixels,
            force_draw_on_top
            );
}
#endif

////////////////////////////////////////////////////////////////////////////////
// #include "easy_soup.cpp"/////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#define ESO_MAX_VERTICES 999999

void _eso_init() {
    COW0._eso_vertex_positions = (real *) calloc(ESO_MAX_VERTICES, 3 * sizeof(real));
    COW0._eso_vertex_colors = (real *) calloc(ESO_MAX_VERTICES, 4 * sizeof(real));
}

void _eso_begin(real *PVM, int primitive, real size_in_pixels, bool use_world_units_instead_of_pixels, bool force_draw_on_top) {
    ASSERT(!COW1._eso_called_eso_begin_before_calling_eso_vertex_or_eso_end);
    COW1._eso_called_eso_begin_before_calling_eso_vertex_or_eso_end = true;
    COW1._eso_primitive = primitive;
    COW1._eso_size_in_pixels = size_in_pixels;
    COW1._eso_use_world_units_instead_of_pixels = use_world_units_instead_of_pixels;
    COW1._eso_overlay = force_draw_on_top;
    COW1._eso_num_vertices = 0;
    ASSERT(PVM);
    memcpy(COW1._eso_PVM, PVM, 16 * sizeof(real));
}

void eso_end() {
    ASSERT(COW1._eso_called_eso_begin_before_calling_eso_vertex_or_eso_end);
    COW1._eso_called_eso_begin_before_calling_eso_vertex_or_eso_end = false;
    _soup_draw(
            COW1._eso_PVM,
            COW1._eso_primitive,
            _SOUP_XYZ,
            _SOUP_RGBA,
            COW1._eso_num_vertices,
            COW0._eso_vertex_positions,
            COW0._eso_vertex_colors,
            0.0,
            0.0,
            0.0,
            0.0,
            COW1._eso_size_in_pixels,
            COW1._eso_use_world_units_instead_of_pixels,
            COW1._eso_overlay
            );
}

void eso_vertex(real x, real y, real z = 0.0) {
    ASSERT(COW1._eso_called_eso_begin_before_calling_eso_vertex_or_eso_end);
    ASSERT(COW1._eso_called_eso_color_at_least_once_before_calling_eso_vertex);
    ASSERT(COW1._eso_num_vertices < ESO_MAX_VERTICES);
    real p[3] = { x, y, z };
    memcpy(COW0._eso_vertex_positions + 3 * COW1._eso_num_vertices, p, 3 * sizeof(real));
    memcpy(COW0._eso_vertex_colors + 4 * COW1._eso_num_vertices, COW1._eso_current_color, 4 * sizeof(real));
    ++COW1._eso_num_vertices;
}

void eso_color(real r, real g, real b, real a = 1.0) {
    COW1._eso_called_eso_color_at_least_once_before_calling_eso_vertex = true;
    COW1._eso_current_color[0] = r;
    COW1._eso_current_color[1] = g;
    COW1._eso_current_color[2] = b;
    COW1._eso_current_color[3] = a;
}

#ifdef SNAIL_CPP
void eso_begin(mat4 PVM, int primitive, real size_in_pixels = 0, bool use_world_units_instead_of_pixels = false, bool force_draw_on_top = false) {
    _eso_begin(PVM.data, primitive, size_in_pixels, use_world_units_instead_of_pixels, force_draw_on_top);
}

void eso_vertex(vec2 xy) {
    eso_vertex(xy[0], xy[1]);
}

void eso_vertex(vec3 xyz) {
    eso_vertex(xyz[0], xyz[1], xyz[2]);
}

void eso_color(vec3 rgb, real a = 1.0) {
    eso_color(rgb[0], rgb[1], rgb[2], a);
}

#endif

////////////////////////////////////////////////////////////////////////////////
// #include "text.cpp"//////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// TODO no_billboard

void _text_draw(
        real *PV,
        char *text,
        real x_world,
        real y_world,
        real z_world,
        real r,
        real g,
        real b,
        real a,
        real font_size_in_pixels,
        real dx_in_pixels,
        real dy_in_pixels,
        bool force_draw_on_top
        ) {
    ASSERT(PV);
    ASSERT(text);

    real window_width_in_pixels, window_height_in_pixels;
    _window_get_size(&window_width_in_pixels, &window_height_in_pixels);

    if (IS_ZERO(font_size_in_pixels)) font_size_in_pixels = 24;

    font_size_in_pixels *= COW0._window_macbook_retina_scale;
    dx_in_pixels *= COW0._window_macbook_retina_scale;
    dy_in_pixels *= COW0._window_macbook_retina_scale;

    static char buffer[99999]; // ~500 chars
    int num_quads = stb_easy_font_print(0, 0, text, NULL, buffer, sizeof(buffer));
    int num_vertices = 4 * num_quads;
    static real vertex_positions[99999];

    char *read_head = buffer;
    for (int k = 0; k < num_vertices; ++k) {
        for (int d = 0; d < 2; ++d) {
            vertex_positions[2 * k + d] = ((float *) read_head)[d];
        }
        read_head += (3 * sizeof(float) + 4);
    }

    real s_NDC[4] = {}; {
        real s_world[4] = { x_world, y_world, z_world, 1 };
        _linalg_mat4_times_vec4_persp_divide(s_NDC, PV, s_world);
    }

    if (IS_BETWEEN(s_NDC[2], -1, 1)) {
        real transform[16] = {}; {
            // Translation(s_NDC) * app_NDC_from_Screen() * Translation(ds_Screen + dims / 2) * Scaling(size, size);

            real TS[16] = {
                font_size_in_pixels / 12, 0, 0, dx_in_pixels + window_width_in_pixels / 2,
                0, font_size_in_pixels / 12, 0, dy_in_pixels + window_height_in_pixels / 2,
                0, 0, 1, 0,
                0, 0, 0, 1,
            };

            _linalg_mat4_times_mat4(transform, (real *) &cow.NDC_from_Screen, TS);
            for (int d = 0; d < 3; ++d) transform[4 * d + 3] += s_NDC[d];
        }

        _soup_draw(
                transform,
                SOUP_QUADS,
                _SOUP_XY,
                _SOUP_RGBA,
                num_vertices,
                vertex_positions,
                NULL,
                r,
                g,
                b,
                a,
                0,
                false,
                force_draw_on_top
                );
    }
}

#ifdef SNAIL_CPP
template<int D_pos = 3, int D_color = 3> void text_draw(
        mat4 PV,
        char *text,
        Vec<D_pos> s_world,
        Vec<D_color> color = { 1.0, 1.0, 1.0 },
        real font_size_in_pixels = 0,
        vec2 ds_in_pixels = { 0.0, 0.0 },
        bool force_draw_on_top = false
        ) {
    STATIC_ASSERT(D_pos == 2 || D_pos == 3 || D_pos == 4);
    STATIC_ASSERT(D_color == 3 || D_color == 4);
    _text_draw(
            PV.data,
            text,
            s_world[0],
            s_world[1],
            D_pos == 3 ? s_world[2] : 0.0,
            color[0],
            color[1],
            color[2],
            D_color == 4 ? color[3] : 1.0,
            font_size_in_pixels,
            ds_in_pixels[0],
            ds_in_pixels[1],
            force_draw_on_top
            );
}
#endif

////////////////////////////////////////////////////////////////////////////////
// #include "gui.cpp"///////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void _gui_begin_frame() {
    COW1._gui_x_curr = 32;
    COW1._gui_y_curr = 48;
    if (cow.key_pressed[config.hotkeys_gui_hide]) {
        COW1._gui_hide_and_disable = !COW1._gui_hide_and_disable;
    }
}

void gui_printf(const char *format, ...) {
    if (COW1._gui_hide_and_disable) { return; }
    static char _text[256] = {};
    {
        va_list arg;
        va_start(arg, format);
        vsnprintf(_text, sizeof(_text), format, arg);
        va_end(arg);
    }

    char *text = _text;
    char *sep = strchr(text, '`'); // fornow hacking in two color text
    if (!sep) {
        _text_draw((real *) &cow.NDC_from_Screen, text, COW1._gui_x_curr, COW1._gui_y_curr, 0.0, 1.0, 1.0, 1.0, 1.0, 0, 0.0, 0.0, true);
    } else {
        real tmp = COW1._gui_x_curr; {
            *sep = 0;
            _text_draw((real *) &cow.NDC_from_Screen, text, COW1._gui_x_curr, COW1._gui_y_curr, 0.0, 1.0, 1.0, 1.0, 1.0, 0, 0.0, 0.0, true);
            COW1._gui_x_curr += 2 * stb_easy_font_width(text);
            text = sep + 1;
            _text_draw((real *) &cow.NDC_from_Screen, text, COW1._gui_x_curr, COW1._gui_y_curr, 0.0, 0.0, 1.0, 1.0, 1.0, 0, 0.0, 0.0, true);
        } COW1._gui_x_curr = tmp;
    }

    COW1._gui_y_curr += 28;
}

char *_gui_hotkey2string(int hotkey) {
    if (hotkey == COW_KEY_TAB) {
        return "TAB";
    }
    if (hotkey == ' ') {
        return "SPACE";
    }
    if (hotkey == COW_KEY_ARROW_LEFT) {
        return "<--";
    }
    if (hotkey == COW_KEY_ARROW_RIGHT) {
        return "-->";
    }
    static char dummy[512];
    if (dummy[2] == 0) {
        for (int i = 0; i < 256; ++i) {
            dummy[2 * i] = char(i);
        }
    }
    return dummy + 2 * hotkey;
}

bool gui_button(char *t, int hotkey = 0) {
    if (COW1._gui_hide_and_disable) { return false; }
    real s_mouse[2];
    _input_get_mouse_position_and_change_in_position_in_world_coordinates((real *) &cow.NDC_from_Screen, s_mouse, s_mouse + 1, NULL, NULL);

    // fornow
    static char text[256];
    if (hotkey) {
        snprintf(text, sizeof(text), "%s `%s", t, _gui_hotkey2string(hotkey));
    } else {
        strcpy(text, t);
    }
    real L = (2 * stb_easy_font_width(text) + 16); // fornow
    if (hotkey) L -= 12;

    real H = 24;
    real box[8] = {
        COW1._gui_x_curr    , COW1._gui_y_curr    ,
        COW1._gui_x_curr + L, COW1._gui_y_curr    ,
        COW1._gui_x_curr + L, COW1._gui_y_curr + H,
        COW1._gui_x_curr    , COW1._gui_y_curr + H,
    };

    bool hot = IS_BETWEEN(s_mouse[0], box[0], box[2]) && IS_BETWEEN(s_mouse[1], box[1], box[5]);

    if (!COW1._gui_selected && ((hot && cow.mouse_left_pressed) || cow.key_pressed[hotkey])) {
        COW1._gui_selected = t;
    }
    if (COW1._gui_selected == t) {
        if (cow.mouse_left_released || cow.key_released[hotkey]) {
            COW1._gui_selected = NULL;
        }
    }

    real r = (COW1._gui_selected != t) ? 0 : .8;
    if (COW1._gui_selected != t) {
        real nudge = SGN(.5 - r) * .1;
        r += nudge; 
        if (hot || cow.key_held[hotkey]) r += nudge; 
    }
    {
        _soup_draw((real *) &cow.NDC_from_Screen, SOUP_QUADS, _SOUP_XY, _SOUP_RGB, 4, box, NULL, r, r, r, 1, 0, false, true);
        _soup_draw((real *) &cow.NDC_from_Screen, SOUP_LINE_LOOP, _SOUP_XY, _SOUP_RGB, 4, box, NULL, 1, 1, 1, 1, 4, false, true);
    }
    COW1._gui_x_curr += 8;
    COW1._gui_y_curr += 4;
    gui_printf(text);
    COW1._gui_y_curr += 8;
    COW1._gui_x_curr -= 8;

    return (COW1._gui_selected == t) && (cow.mouse_left_pressed || cow.key_pressed[hotkey]);
}

void gui_checkbox(char *name, bool *t, int hotkey = 0) {
    if (COW1._gui_hide_and_disable) { return; }
    real s_mouse[2];
    _input_get_mouse_position_and_change_in_position_in_world_coordinates((real *) &cow.NDC_from_Screen, s_mouse, s_mouse + 1, NULL, NULL);
    real L = 16;
    real box[8] = {
        COW1._gui_x_curr    , COW1._gui_y_curr    ,
        COW1._gui_x_curr + L, COW1._gui_y_curr    ,
        COW1._gui_x_curr + L, COW1._gui_y_curr + L,
        COW1._gui_x_curr    , COW1._gui_y_curr + L,
    };
    bool hot = IS_BETWEEN(s_mouse[0], box[0], box[2]) && IS_BETWEEN(s_mouse[1], box[1], box[5]);

    if ((hot && cow.mouse_left_pressed) || cow.key_pressed[hotkey]) {
        *t = !(*t);
        if (!COW1._gui_selected) COW1._gui_selected = t;
    }
    if (COW1._gui_selected == t) {
        if (cow.mouse_left_released || cow.key_released[hotkey]) {
            COW1._gui_selected = NULL;
        }
    }

    real r = (!*t) ? 0 : 1;
    if (COW1._gui_selected != t) {
        real nudge = SGN(.5 - r) * .1;
        r += nudge; 
        if (hot || cow.key_held[hotkey]) r += nudge; 
    }
    {
        _soup_draw((real *) &cow.NDC_from_Screen, SOUP_QUADS, _SOUP_XY, _SOUP_RGB, 4, box, NULL, r, r, r, 1, 0, false, true);
        _soup_draw((real *) &cow.NDC_from_Screen, SOUP_LINE_LOOP, _SOUP_XY, _SOUP_RGB, 4, box, NULL, 1, 1, 1, 1, 4, false, true);
    }
    COW1._gui_x_curr += 2 * L;
    if (hotkey) {
        gui_printf("%s `%s", name, _gui_hotkey2string(hotkey));
    } else {
        gui_printf(name);
    }
    COW1._gui_x_curr -= 2 * L;
}

void gui_readout(char *name, int *t) {
    if (!name) name = "";
    char *join = (char *)((name) ? " " : "");
    gui_printf("%s%s%d", name, join, *t);
}

void gui_readout(char *name, real *t) {
    if (!name) name = "";
    char *join = (char *)((name) ? " " : "");
    gui_printf("%s%s%.2lf", name, join, *t);
}

void _gui_slider(char *text, void *t, real *t_copy, real a, real b) {
    COW1._gui_y_curr += 8;
    real s_mouse[2];
    _input_get_mouse_position_and_change_in_position_in_world_coordinates((real *) &cow.NDC_from_Screen, s_mouse, s_mouse + 1, NULL, NULL);
    real w = 166;
    real band[4] = { COW1._gui_x_curr, COW1._gui_y_curr, COW1._gui_x_curr + w, COW1._gui_y_curr };
    real s_dot[2] = { LERP(INVERSE_LERP(*t_copy, a, b), band[0], band[2]), band[1] };

    if (cow._mouse_owner == COW_MOUSE_OWNER_NONE || cow._mouse_owner == COW_MOUSE_OWNER_GUI) {
        bool is_near = _linalg_vecX_squared_distance(2, s_dot, s_mouse) < COW0._window_macbook_retina_scale * 16;
        if (is_near) {
            COW1._gui_hot = t;
            cow._mouse_owner = COW_MOUSE_OWNER_GUI;
        }
        if (COW1._gui_hot == t && !is_near) {
            COW1._gui_hot = 0;
            if (COW1._gui_selected != t) cow._mouse_owner = COW_MOUSE_OWNER_NONE;
        }
    }
    if (!COW1._gui_selected && (COW1._gui_hot == t) && cow.mouse_left_pressed) {
        cow._mouse_owner = COW_MOUSE_OWNER_GUI;
        COW1._gui_selected = t;
    }
    if (COW1._gui_selected == t) {
        if (cow.mouse_left_held) {
            *t_copy = LERP(CLAMP(INVERSE_LERP(s_mouse[0], band[0], band[2]), 0, 1), a, b);
        }
        if (cow.mouse_left_released) {
            if (COW1._gui_hot != t) cow._mouse_owner = COW_MOUSE_OWNER_NONE;
            COW1._gui_selected = 0;
        }
    }
    {
        _soup_draw((real *) &cow.NDC_from_Screen, SOUP_LINES, _SOUP_XY, _SOUP_RGB, 2, band, NULL, .6, .6, .6, 1, 6, false, true);
        real r = (COW1._gui_selected == t) ? 1 : (COW1._gui_hot == t) ? .9 : .8;
        _soup_draw((real *) &cow.NDC_from_Screen, SOUP_POINTS, _SOUP_XY, _SOUP_RGB, 1, s_dot, NULL, r, r, r, 1, (COW1._gui_hot == t && COW1._gui_selected != t) ? 17 : 14, false, true);
    }
    COW1._gui_y_curr -= 8;
    COW1._gui_x_curr += w + 16;
    gui_printf(text);
    COW1._gui_x_curr -= w + 16;
}

void gui_slider(char *name, int *t, int a, int b, int j, int k, bool loop) {
    if (COW1._gui_hide_and_disable) { return; }
    real tmp = real(*t);
    static char text[256]; {
        if (!j && !k) {
            snprintf(text, sizeof(text), "%s", name);
        } else {
            snprintf(text, sizeof(text), "%s %d `%s %s", name, *t, j ? _gui_hotkey2string(j) : "~", k ? _gui_hotkey2string(k) : "~");
        }
    }
    _gui_slider(text, t, &tmp, a, b);
    *t = int(.5 + tmp);
    if (cow.key_pressed[k]) ++(*t);
    if (cow.key_pressed[j]) --(*t);
    if (cow.key_pressed[k] || cow.key_pressed[j]) {
        *t = (!loop) ? CLAMP(*t, a, b) : a + MODULO(*t - a, (b + 1) - a);
    }
}

void gui_slider(char *name, real *t, real a, real b, bool readout_in_degrees) {
    if (COW1._gui_hide_and_disable) { return; }
    static char text[256];
    if (!readout_in_degrees) {
        snprintf(text, sizeof(text), "%s %.1lf", name, *t);
    } else {
        int tmp = (int) round(DEG(*((real *) t)));
        snprintf(text, sizeof(text), "%s %d deg", name, tmp);
    }

    _gui_slider(text, t, t, a, b);
}

////////////////////////////////////////////////////////////////////////////////
// #include "indexed.cpp"///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void _itri_init() {
    COW0._itri_shader_program = shader_compile_and_build_program(COWX._itri_vert, COWX._itri_frag);
    glGenVertexArrays(1, &COW0._itri_VAO);
    glGenBuffers(NELEMS(COW0._itri_VBO), COW0._itri_VBO);
    glGenBuffers(1, &COW0._itri_EBO);
    stbi_set_flip_vertically_on_load(true);
}

bool _itri_texture_find(int *i_texture, char *texture_filename) { // fornow O(n)
    bool already_loaded = false;
    for ((*i_texture) = 0; (*i_texture) < COW1._itri_num_textures; (*i_texture)++) {
        if (strcmp(texture_filename, COW1._itri_texture_filenames[*i_texture]) == 0) {
            already_loaded = true;
            break;
        }
    }
    return already_loaded;
}

int _itri_texture_get_format(int nrChannels) {
    ASSERT(nrChannels == 1 || nrChannels == 3 || nrChannels == 4);
    return (nrChannels == 1) ? GL_RED : (nrChannels == 3) ? GL_RGB : GL_RGBA;
}

int _itri_texture_create(char *texture_filename, int width, int height, int nrChannels, u8 *data) {
    ASSERT(texture_filename);
    ASSERT(COW1._itri_num_textures < ITRI_MAX_NUM_TEXTURES);
    ASSERT(strlen(texture_filename) < ITRI_MAX_FILENAME_LENGTH);
    ASSERT(width > 0);
    ASSERT(height > 0);
    ASSERT(data);

    int format = _itri_texture_get_format(nrChannels);

    int i_texture = COW1._itri_num_textures++;
    strcpy(COW1._itri_texture_filenames[i_texture], texture_filename);

    glGenTextures(1, COW1._itri_textures + i_texture);
    glActiveTexture(GL_TEXTURE0 + i_texture);
    glBindTexture(GL_TEXTURE_2D, COW1._itri_textures[i_texture]);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);   
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, format, GL_UNSIGNED_BYTE, data);
    glGenerateMipmap(GL_TEXTURE_2D);

    return i_texture;
}

void _itri_texture_update(char *texture_filename, int width, int height, int nrChannels, u8 *data) {
    ASSERT(texture_filename);
    ASSERT(width > 0);
    ASSERT(height > 0);
    ASSERT(data);

    int format = _itri_texture_get_format(nrChannels);

    int i_texture;
    ASSERT(_itri_texture_find(&i_texture, texture_filename));
    glActiveTexture(GL_TEXTURE0 + i_texture);
    glBindTexture(GL_TEXTURE_2D, COW1._itri_textures[i_texture]);

    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, format, GL_UNSIGNED_BYTE, data);
}

int _itri_texture_load(char *texture_filename) {
    ASSERT(texture_filename);
    int width, height, nrChannels;
    u8 *data = stbi_load(texture_filename, &width, &height, &nrChannels, 0);
    ASSERT(data);
    ASSERT(width > 0);
    ASSERT(height > 0);
    ASSERT(nrChannels == 3 || nrChannels == 4);
    int i_texture = _itri_texture_create(texture_filename, width, height, nrChannels, data);
    stbi_image_free(data);
    return i_texture;
}

struct Texture {
    char *filename;
    int width;
    int height;
    int nrChannels;
    u8 *data;
};

int itri_texture_create(Texture *texture) {
    return _itri_texture_create(texture->filename, texture->width, texture->height, texture->nrChannels, texture->data);
}

void itri_texture_update(Texture *texture) {
    _itri_texture_update(texture->filename, texture->width, texture->height, texture->nrChannels, texture->data);
}

void _itri_draw(
        real *P,
        real *V,
        real *M,
        int num_triangles,
        int *triangle_indices,
        int num_vertices,
        real *vertex_positions,
        real *vertex_normals = NULL,
        real *vertex_colors = NULL,
        real r_if_vertex_colors_is_NULL = 1,
        real g_if_vertex_colors_is_NULL = 1,
        real b_if_vertex_colors_is_NULL = 1,
        real *vertex_texCoords = NULL,
        char *texture_filename = NULL
        ) {
    if (num_triangles == 0) { return; } // NOTE: num_triangles zero is now valid input
    ASSERT(P);
    ASSERT(V);
    ASSERT(M);
    ASSERT(vertex_positions);
    real color_if_vertex_colors_is_NULL[3] = { r_if_vertex_colors_is_NULL, g_if_vertex_colors_is_NULL, b_if_vertex_colors_is_NULL };

    int i_texture = -1;
    if (texture_filename) {
        if (!_itri_texture_find(&i_texture, texture_filename)) {
            i_texture = _itri_texture_load(texture_filename);
        }
    }

    glBindVertexArray(COW0._itri_VAO);
    int i_attrib = 0;
    auto guarded_push = [&](void *array, int dim) {
        glDisableVertexAttribArray(i_attrib); // fornow
        if (array) {
            glBindBuffer(GL_ARRAY_BUFFER, COW0._itri_VBO[i_attrib]);
            glBufferData(GL_ARRAY_BUFFER, num_vertices * dim * sizeof(real), array, GL_DYNAMIC_DRAW);
            glVertexAttribPointer(i_attrib, dim, GL_REAL, 0, 0, NULL);
            glEnableVertexAttribArray(i_attrib);
        }
        ++i_attrib;
    };

    guarded_push(vertex_positions, 3);
    guarded_push(vertex_normals, 3);
    guarded_push(vertex_colors, 3);
    guarded_push(vertex_texCoords, 2);

    ASSERT(COW0._itri_shader_program);
    glUseProgram(COW0._itri_shader_program);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, COW0._itri_EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3 * num_triangles * sizeof(u32), triangle_indices, GL_DYNAMIC_DRAW);

    _shader_set_uniform_mat4(COW0._itri_shader_program, "P", P);
    _shader_set_uniform_mat4(COW0._itri_shader_program, "V", V);
    _shader_set_uniform_mat4(COW0._itri_shader_program, "M", M);
    { // fornow scavenge the camera position from V
        real C[16];
        _linalg_mat4_inverse(C, V);
        real eye_World[4] = { C[3], C[7], C[11], 1 };
        _shader_set_uniform_vec4(COW0._itri_shader_program, "eye_World", eye_World);
    }
    _shader_set_uniform_bool(COW0._itri_shader_program, "has_vertex_colors", vertex_colors != NULL);
    _shader_set_uniform_bool(COW0._itri_shader_program, "has_vertex_normals", vertex_normals != NULL);
    _shader_set_uniform_bool(COW0._itri_shader_program, "has_vertex_texCoords", vertex_texCoords != NULL);
    _shader_set_uniform_bool(COW0._itri_shader_program, "has_texture", texture_filename != NULL);
    _shader_set_uniform_int (COW0._itri_shader_program, "i_texture", MAX(0, i_texture));
    _shader_set_uniform_vec3(COW0._itri_shader_program, "color_if_vertex_colors_is_NULL", color_if_vertex_colors_is_NULL);

    glDrawElements(GL_TRIANGLES, 3 * num_triangles, GL_UNSIGNED_INT, NULL);
}

#ifdef SNAIL_CPP
void itri_draw(
        mat4 P,
        mat4 V,
        mat4 M,
        int num_triangles,
        int3 *triangle_indices,
        int num_vertices,
        vec3 *vertex_positions,
        vec3 *vertex_normals,
        vec3 *vertex_colors,
        vec3 color_if_vertex_colors_is_NULL,
        vec2 *vertex_texCoords,
        char *texture_filename
        ) {
    _itri_draw(
            P.data,
            V.data,
            M.data,
            num_triangles,
            (int *) triangle_indices,
            num_vertices,
            (real *) vertex_positions,
            (real *) vertex_normals,
            (real *) vertex_colors,
            color_if_vertex_colors_is_NULL[0],
            color_if_vertex_colors_is_NULL[1],
            color_if_vertex_colors_is_NULL[2],
            (real *) vertex_texCoords,
            texture_filename
            );
}
#endif

////////////////////////////////////////////////////////////////////////////////
// #include "camera.cpp"////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// this is an "ortho camera" that points at (o_x, o_y)   
struct Camera2D {
    real screen_height_World;
    real o_x;
    real o_y;
};

// this is an "orbit camera" that points at the origin*  
//                               *unless o_x, o_y nonzero
//                                                       
// it is stuck to the surface of a sphere at (theta, phi)
// these angles can also be interpreted as (yaw, pitch)  
//                                                       
// the sphere radius (camera distance) is given by*      
//     (screen_height_World / 2) / tan(angle_of_view / 2)
// *we define it this way for compatability with Camera2D
struct Camera3D {
    union {
        real persp_distance_to_origin;
        real ortho_screen_height_World;
    };
    real angle_of_view; // set to 0 to use an ortho camera
    real theta; // yaw
    real phi; // pitch
    real o_x;
    real o_y;
};

void camera_attach_to_gui(Camera2D *camera) {
    gui_slider("screen_height_World", &camera->screen_height_World, 0.1, 100.0, false);
    gui_slider("o_x", &camera->o_x, -camera->screen_height_World, camera->screen_height_World, false);
    gui_slider("o_y", &camera->o_y, -camera->screen_height_World, camera->screen_height_World, false);
}

void camera_attach_to_gui(Camera3D *camera) {
    if (IS_ZERO(camera->angle_of_view)) { // ortho
        gui_slider("ortho_screen_height_World", &camera->ortho_screen_height_World, 0.1, 100.0, false);
    } else {
        gui_slider("persp_distance_to_origin", &camera->persp_distance_to_origin, 0.1, 100.0, false);
    }
    gui_slider("angle_of_view", &camera->angle_of_view, RAD(0.0), RAD(180.0), true);
    gui_slider("theta", &camera->theta, RAD(-180.0), RAD(180.0), true);
    gui_slider("phi", &camera->phi, RAD(-90.0), RAD(90.0), true);
    gui_slider("o_x", &camera->o_x, -camera->persp_distance_to_origin, camera->persp_distance_to_origin, false);
    gui_slider("o_y", &camera->o_y, -camera->persp_distance_to_origin, camera->persp_distance_to_origin, false);
}

void _camera_get_PV(Camera2D *camera, real *PV) {
    _window_get_PV_ortho(PV, camera->screen_height_World, camera->o_x, camera->o_y);
}

real camera_get_screen_height_World(Camera3D *camera) {
    if (IS_ZERO(camera->angle_of_view)) { // ortho
        return camera->ortho_screen_height_World;
    } else { // persp
        real theta = camera->angle_of_view / 2;
        return 2 * tan(theta) * camera->persp_distance_to_origin;
    }
}

void _camera_get_coordinate_system(Camera3D *camera, real *C_out, real *x_out = NULL, real *y_out = NULL, real *z_out = NULL, real *o_out = NULL, real *R_out = NULL) {
    real camera_o_z = camera->persp_distance_to_origin;

    real C[16]; {
        real T[16] = {
            1, 0, 0, camera->o_x,
            0, 1, 0, camera->o_y, // fornow
            0, 0, 1, camera_o_z,
            0, 0, 0, 1,
        };
        real c_x = cos(camera->phi);
        real s_x = sin(camera->phi);
        real c_y = cos(camera->theta);
        real s_y = sin(camera->theta);
        real R_y[16] = {
            c_y , 0, s_y, 0,
            0   , 1, 0  , 0,
            -s_y, 0, c_y, 0,
            0   , 0, 0  , 1,
        };
        real R_x[16] = {
            1,   0,    0, 0, 
            0, c_x, -s_x, 0,
            0, s_x,  c_x, 0,
            0  , 0   , 0, 1,
        };
        _linalg_mat4_times_mat4(C,      R_x, T);
        _linalg_mat4_times_mat4(C, R_y, C     );
    }

    real xyzo[4][3] = {}; {
        for (int c = 0; c < 4; ++c) for (int r = 0; r < 3; ++r) xyzo[c][r] = _LINALG_4X4(C, r, c);
    }
    if (C_out) memcpy(C_out, C, sizeof(C));
    {
        if (x_out) memcpy(x_out, xyzo[0], 3 * sizeof(real));
        if (y_out) memcpy(y_out, xyzo[1], 3 * sizeof(real));
        if (z_out) memcpy(z_out, xyzo[2], 3 * sizeof(real));
        if (o_out) memcpy(o_out, xyzo[3], 3 * sizeof(real));
    }
    if (R_out) {
        real tmp[] = {
            C[0], C[1], C[ 2], 0,
            C[4], C[5], C[ 6], 0,
            C[8], C[9], C[10], 0,
            0   , 0   , 0    , 1,
        };
        memcpy(R_out, tmp, sizeof(tmp));
    }
}

void _camera_get_P(Camera3D *camera, real *P) {
    if (IS_ZERO(camera->angle_of_view)) {
        _window_get_P_ortho(P, camera->ortho_screen_height_World);
    } else {
        _window_get_P_perspective(P, camera->angle_of_view);
    }
}

void _camera_get_V(Camera3D *camera, real *V) {
    _camera_get_coordinate_system(camera, V);
    _linalg_mat4_inverse(V, V);
}

void _camera_get_PV(Camera3D *camera, real *PV) {
    ASSERT(PV);
    real P[16], V[16];
    _camera_get_P(camera, P);
    _camera_get_V(camera, V);
    _linalg_mat4_times_mat4(PV, P, V);
}

void camera_move(Camera2D *camera, bool disable_pan = false, bool disable_zoom = false) {
    real NDC_from_World[16] = {};
    _camera_get_PV(camera, NDC_from_World);
    real x, y, dx, dy;
    _input_get_mouse_position_and_change_in_position_in_world_coordinates(NDC_from_World, &x, &y, &dx, &dy);
    if (!disable_pan && cow.mouse_right_held) {
        camera->o_x -= dx;
        camera->o_y -= dy;
    }
    else if (!disable_zoom && !IS_ZERO(cow.mouse_wheel_offset)) {
        camera->screen_height_World *= (1 - .1 * cow.mouse_wheel_offset);

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

        real mouse_World[4] = { x, y, 0, 1 };
        real a[4] = {};
        real inv_P_prime[16] = {}; {
            Camera2D tmp = { camera->screen_height_World };
            _camera_get_PV(&tmp, inv_P_prime);
            _linalg_mat4_inverse(inv_P_prime, inv_P_prime);
        }
        real mouse_NDC[4] = {}; {
            _linalg_mat4_times_vec4_persp_divide(mouse_NDC, NDC_from_World, mouse_World); 
        }
        _linalg_mat4_times_vec4_persp_divide(a, inv_P_prime, mouse_NDC);
        for (int d = 0; d < 2; ++d) (&camera->o_x)[d] = mouse_World[d] - a[d];
    }
}

void camera_move(Camera3D *camera, bool disable_pan = false, bool disable_zoom = false, bool disable_rotate = false) {
    disable_rotate |= (cow._mouse_owner != COW_MOUSE_OWNER_NONE);
    { // 2D transforms
        Camera2D tmp2D = { camera_get_screen_height_World(camera), camera->o_x, camera->o_y };
        camera_move(&tmp2D, disable_pan, disable_zoom);

        if (IS_ZERO(camera->angle_of_view)) { // ortho
            camera->ortho_screen_height_World = tmp2D.screen_height_World;
        } else { // persp
            //                   r_y
            //               -   |  
            //            -      |  
            //         -         |  
            //      - ) theta    |  
            // eye-------------->o  
            //         D            
            real r_y = tmp2D.screen_height_World / 2;
            real theta = camera->angle_of_view / 2;
            camera->persp_distance_to_origin = r_y / tan(theta);
        }

        camera->o_x = tmp2D.o_x;
        camera->o_y = tmp2D.o_y;
    }
    if (!disable_rotate && cow.mouse_left_held) {
        real fac = 2;
        camera->theta -= fac * cow.mouse_change_in_position_NDC[0];
        camera->phi += fac * cow.mouse_change_in_position_NDC[1];
        camera->phi = CLAMP(camera->phi, -RAD(90), RAD(90));
    }
}

#ifdef SNAIL_CPP
struct OrthogonalCoordinateSystem3D {
    //   = [| | | |]   [\ | / |]
    // C = [x y z o] = [- R - o]
    //   = [| | | |]   [/ | \ |]
    //   = [0 0 0 1]   [- 0 - 1]
    mat4 C;
    vec3 x;
    vec3 y;
    vec3 z;
    vec3 o;
    mat4 R; // stored as [\ | / |]
    //           [- R - 0]
    //           [/ | \ |]
    //           [- 0 - 1]
};

mat4 camera_get_PV(Camera2D *camera) {
    mat4 ret;
    _camera_get_PV(camera, ret.data);
    return ret;
}

OrthogonalCoordinateSystem3D camera_get_coordinate_system(Camera3D *camera) {
    mat4 C;
    vec3 x, y, z, o;
    mat4 R;
    _camera_get_coordinate_system(camera, (real *) &C, (real *) &x, (real *) &y, (real *) &z, (real *) &o, (real *) &R);
    return { C, x, y, z, o, R };
}

vec3 camera_get_origin(Camera3D *camera) {
    vec3 o;
    _camera_get_coordinate_system(camera, NULL, NULL, NULL, NULL, (real *) &o, NULL);
    return o;
}

mat4 camera_get_P(Camera3D *camera) {
    mat4 ret;
    _camera_get_P(camera, ret.data);
    return ret;
}

mat4 camera_get_V(Camera3D *camera) {
    mat4 ret;
    _camera_get_V(camera, ret.data);
    return ret;
}

mat4 camera_get_C(Camera3D *camera) {
    mat4 ret;
    _camera_get_coordinate_system(camera, ret.data);
    return ret;
}

mat4 camera_get_PV(Camera3D *camera) {
    mat4 ret;
    _camera_get_PV(camera, ret.data);
    return ret;
}
#endif

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// extras //////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// #include "color.cpp"/////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#ifdef SNAIL_CPP
struct {
    vec3 red    = { 249./255,  38./255, 114./255 }; /** red */
    vec3 orange = { 253./255, 151./255,  31./255 }; /** orange */
    vec3 yellow = { 255./255, 255./255,  50./255 }; /** not the actual monokai yellow cause i don't like it */
    vec3 green  = { 166./255, 226./255,  46./255 }; /** green */
    vec3 blue   = { 102./255, 217./255, 239./255 }; /** blue */
    vec3 purple = { 174./255, 129./255, 255./255 }; /** purple */
    vec3 white  = { 255./255, 255./255, 255./255 }; /** full white */
    vec3 gray   = { 127./255, 127./255, 127./255 }; /** 50% gray */
    vec3 black  = {   0./255,   0./255,   0./255 }; /** full black */
    vec3 brown  = { 123./255,  63./255,   0./255 }; /** monokai doesn't actually define a brown so i made this one up */
} monokai;

vec3 color_kelly(int i) {
    static vec3 _kelly_colors[]={{255./255,179./255,0./255},{128./255,62./255,117./255},{255./255,104./255,0./255},{166./255,189./255,215./255},{193./255,0./255,32./255},{206./255,162./255,98./255},{129./255,112./255,102./255},{0./255,125./255,52./255},{246./255,118./255,142./255},{0./255,83./255,138./255},{255./255,122./255,92./255},{83./255,55./255,122./255},{255./255,142./255,0./255},{179./255,40./255,81./255},{244./255,200./255,0./255},{127./255,24./255,13./255},{147./255,170./255,0./255},{89./255,51./255,21./255},{241./255,58./255,19./255},{35./255,44./255,22./255}};
    return _kelly_colors[MODULO(i, NELEMS(_kelly_colors))];
}

vec3 color_plasma(real t) {
    t = CLAMP(t, 0.0, 1.0);
    const vec3 c0 = { 0.05873234392399702, 0.02333670892565664, 0.5433401826748754 };
    const vec3 c1 = { 2.176514634195958, 0.2383834171260182, 0.7539604599784036 };
    const vec3 c2 = { -2.689460476458034, -7.455851135738909, 3.110799939717086 };
    const vec3 c3 = { 6.130348345893603, 42.3461881477227, -28.51885465332158 };
    const vec3 c4 = { -11.10743619062271, -82.66631109428045, 60.13984767418263 };
    const vec3 c5 = { 10.02306557647065, 71.41361770095349, -54.07218655560067 };
    const vec3 c6 = { -3.658713842777788, -22.93153465461149, 18.19190778539828 };
    return c0+t*(c1+t*(c2+t*(c3+t*(c4+t*(c5+t*c6)))));
}

vec3 color_rainbow_swirl(real t) {
    #define Q(o) (.5 + .5 * cos(6.28 * ((o) - t)))
    return { Q(0), Q(.33), Q(-.33) };
    #undef Q
}
#endif

////////////////////////////////////////////////////////////////////////////////
// #include "sbuff.cpp"/////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <typename T> struct StretchyBuffer {
    int length;
    int capacity;
    T *data;
    T &operator [](int index) { return data[index]; }
};

template <typename T> void sbuff_push_back(StretchyBuffer<T> *buffer, T element) {
    if (buffer->capacity == 0) {
        buffer->capacity = 16;
        buffer->data = (T *) malloc(buffer->capacity * sizeof(T));
    }
    if (buffer->length == buffer->capacity) {
        buffer->capacity *= 2;
        buffer->data = (T *) realloc(buffer->data, buffer->capacity * sizeof(T));
    }
    buffer->data[buffer->length++] = element;
}

template <typename T> void sbuff_free(StretchyBuffer<T> *buffer) {
    buffer->length = 0;
    buffer->capacity = 0;
    if (buffer->data) {
        free(buffer->data);
    }
    buffer->data = NULL;
}

////////////////////////////////////////////////////////////////////////////////
// #include "mesh.cpp"//////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#ifdef SNAIL_CPP
struct Soup3D {
    int primitive;
    int dimension_of_positions;
    int dimension_of_colors;
    int num_vertices;
    vec3 *vertex_positions;
    vec3 *vertex_colors;

    void draw(
            mat4 PVM,
            vec3 color_if_vertex_colors_is_NULL,
            real size_in_pixels,
            bool force_draw_on_top);

    void _dump_for_meshlib(char *filename, char *name);
};

struct IndexedTriangleMesh3D {
    int num_vertices;
    int num_triangles;
    vec3 *vertex_positions;
    vec3 *vertex_normals;
    vec3 *vertex_colors;
    int3 *triangle_indices;
    vec2 *vertex_texCoords;
    char *texture_filename;

    void draw(
            mat4 P,
            mat4 V,
            mat4 M,
            vec3 color_if_vertex_colors_is_NULL,
            char *texture_filename_if_texture_filename_is_NULL
            );

    void _dump_for_meshlib(char *filename, char *name);
};

void Soup3D::draw(
        mat4 PVM,
        vec3 color_if_vertex_colors_is_NULL = { 1.0, 0.0, 1.0 },
        real size_in_pixels = 0,
        bool force_draw_on_top = false
        ) {
    soup_draw(
            PVM,
            primitive,
            num_vertices,
            vertex_positions,
            vertex_colors,
            color_if_vertex_colors_is_NULL,
            size_in_pixels,
            force_draw_on_top);
}

void IndexedTriangleMesh3D::draw(
        mat4 P,
        mat4 V,
        mat4 M,
        vec3 color_if_vertex_colors_is_NULL = { 1.0, 0.0, 1.0 },
        char *texture_filename_if_texture_filename_is_NULL = NULL
        ) {
    itri_draw(
            P,
            V,
            M,
            num_triangles,
            triangle_indices,
            num_vertices,
            vertex_positions,
            vertex_normals,
            vertex_colors,
            color_if_vertex_colors_is_NULL,
            vertex_texCoords,
            (texture_filename) ? texture_filename : texture_filename_if_texture_filename_is_NULL
            );
}

void Soup3D::_dump_for_meshlib(char *filename, char *name) {
    ASSERT(primitive == SOUP_TRIANGLE_MESH);
    FILE *fp = fopen(filename, "w");
    fprintf(fp, "const int _meshlib_soup_%s_num_vertices = %d;\n", name, num_vertices); 
    fprintf(fp, "const vec3 _meshlib_soup_%s_vertex_positions[_meshlib_soup_%s_num_vertices] = {\n    ", name, name); for (int i = 0; i < num_vertices; ++i) fprintf(fp, "{%.3lf,%.3lf,%.3lf},",vertex_positions[i][0],vertex_positions[i][1],vertex_positions[i][2]); fprintf(fp, "};\n");
    fclose(fp);
}

void IndexedTriangleMesh3D::_dump_for_meshlib(char *filename, char *name) {
    FILE *fp = fopen(filename, "w");
    fprintf(fp, "const int _meshlib_itri_%s_num_triangles = %d;\n", name, num_triangles); 
    fprintf(fp, "const int _meshlib_itri_%s_num_vertices = %d;\n", name, num_vertices); 
    fprintf(fp, "const int3 _meshlib_itri_%s_triangle_indices[_meshlib_itri_%s_num_triangles] = {\n    ", name, name); for (int i = 0; i < num_triangles; ++i) fprintf(fp, "{%d,%d,%d},",triangle_indices[i].i,triangle_indices[i].j,triangle_indices[i].k); fprintf(fp, "};\n");
    fprintf(fp, "const vec3 _meshlib_itri_%s_vertex_positions[_meshlib_itri_%s_num_vertices] = {\n    ", name, name); for (int i = 0; i < num_vertices; ++i) fprintf(fp, "{%.3lf,%.3lf,%.3lf},",vertex_positions[i][0],vertex_positions[i][1],vertex_positions[i][2]); fprintf(fp, "};\n");
    fprintf(fp, "const vec3 _meshlib_itri_%s_vertex_normals  [_meshlib_itri_%s_num_vertices] = {\n    ", name, name); for (int i = 0; i < num_vertices; ++i) fprintf(fp, "{%.3lf,%.3lf,%.3lf},",vertex_normals[i][0],vertex_normals[i][1],vertex_normals[i][2]); fprintf(fp, "};\n");
    fclose(fp);
}

////////////////////////////////////////////////////////////////////////////////
// #include "_mesh_utils.cpp"///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#if 0
IndexedTriangleMesh3D bunny = _mesh_utils_indexed_triangle_mesh_load("data_fancy_bunny", true, true, false);
bunny._dump_for_meshlib("out.txt", "bunny");
#endif

#if 0
Soup3D bunny = _mesh_utils_soup_TRIANGLES_load("data_basic_bunny", true);
bunny._dump_for_meshlib("out.txt", "bunny");
#endif

void _mesh_utils_transform_vertex_positions_to_double_unit_box(int num_vertices, vec3 *vertex_positions) {
    vec3 L = V3(HUGE, HUGE, HUGE);
    vec3 R = V3(-HUGE, -HUGE, -HUGE);
    for (int i = 0; i < num_vertices; ++i) {
        L = cwiseMin(L, vertex_positions[i]);
        R = cwiseMax(R, vertex_positions[i]);
    }
    vec3 center = .5 * (L + R);
    vec3 size = R - L;
    real largest = MAX(MAX(size[0], size[1]), size[2]);
    for (int i = 0; i < num_vertices; ++i) {
        vertex_positions[i] -= center;
        vertex_positions[i] *= (2 / largest);
    }
}

void _mesh_utils_indexed_triangle_mesh_alloc_compute_and_store_area_weighted_vertex_normals(IndexedTriangleMesh3D *fancy_mesh) {
    ASSERT(fancy_mesh->vertex_normals == NULL);
    if (1) { // () _itri_triangle_mesh_alloc_compute_and_store_area_weighted_vertex_normals
        // TODO allocate fancy_mesh->vertex_normals        
        // TODO write entries of fancy_mesh->vertex_normals
        fancy_mesh->vertex_normals = (vec3 *) calloc(fancy_mesh->num_vertices, sizeof(vec3));
        for (int i_triangle = 0; i_triangle < fancy_mesh->num_triangles; ++i_triangle) {
            int3 ijk = fancy_mesh->triangle_indices[i_triangle];
            real A;
            vec3 n_hat;
            {
                vec3 abc[3];
                for (int d = 0; d < 3; ++d) {
                    abc[d] = fancy_mesh->vertex_positions[ijk[d]];
                }
                vec3 n = cross(abc[1] - abc[0], abc[2] - abc[0]);
                real mag_n = norm(n);
                A = norm(n) / 2;
                n_hat = n / mag_n;
            }
            for (int d = 0; d < 3; ++d) {
                fancy_mesh->vertex_normals[ijk[d]] += A * n_hat;
            }
        }
        for (int i_vertex = 0; i_vertex < fancy_mesh->num_vertices; ++i_vertex) {
            fancy_mesh->vertex_normals[i_vertex] = normalized(fancy_mesh->vertex_normals[i_vertex]);
        }
    }
}

void _mesh_utils_indexed_triangle_mesh_merge_duplicated_vertices(IndexedTriangleMesh3D *fancy_mesh) {
    int new_num_vertices = 0;
    vec3 *new_vertex_positions = (vec3 *) calloc(fancy_mesh->num_vertices, sizeof(vec3)); // (more space than we'll need)
    if (1) { // [] _itri_mesh_merge_duplicated_vertices
        // TODO set new_num_vertices and entries of new_vertex_positions                   
        // TODO overwrite entries of fancy_mesh->triangle_indices with new triangle indices
        // NOTE it is OK if your implementation is slow (mine takes ~5 seconds to run)     
        // NOTE please don't worry about space efficiency at all                           
        int *primal  = (int *) calloc(fancy_mesh->num_vertices, sizeof(int));
        int *new_index = (int *) calloc(fancy_mesh->num_vertices, sizeof(int));
        for (int i = 0; i < fancy_mesh->num_vertices; ++i) {
            primal[i] = -1;
            new_index[i] = -1;
        }
        for (int i = 0; i < fancy_mesh->num_vertices; ++i) {
            if (primal[i] != -1) {
                continue;
            }
            vec3 p_i = fancy_mesh->vertex_positions[i];
            for (int j = i + 1; j < fancy_mesh->num_vertices; ++j) {
                if (primal[j] != -1) {
                    continue;
                }
                vec3 p_j = fancy_mesh->vertex_positions[j];
                if (IS_ZERO(squaredNorm(p_i - p_j))) {
                    ASSERT(primal[j] == -1);
                    primal[j] = i;
                }
            }
        }
        int k = 0;
        for (int i = 0; i < fancy_mesh->num_vertices; ++i) {
            if (primal[i] == -1) {
                ++new_num_vertices;
                new_vertex_positions[k] = fancy_mesh->vertex_positions[i];
                new_index[i] = k;
                ++k;
            }
        }
        for (int i_triangle = 0; i_triangle < fancy_mesh->num_triangles; ++i_triangle) {
            for (int d = 0; d < 3; ++d) {
                int i = fancy_mesh->triangle_indices[i_triangle][d];
                fancy_mesh->triangle_indices[i_triangle][d] = (primal[i] == -1) ? new_index[i] : new_index[primal[i]];
            }
        }
        free(primal);
        free(new_index);
    }
    if (new_num_vertices) {
        fancy_mesh->num_vertices = new_num_vertices;
        free(fancy_mesh->vertex_positions);
        fancy_mesh->vertex_positions = new_vertex_positions;
    }
}

IndexedTriangleMesh3D _mesh_utils_indexed_triangle_mesh_load(char *filename, bool transform_vertex_positions_to_double_unit_box, bool compute_normals, bool merge_duplicated_vertices) {
    IndexedTriangleMesh3D fancy_mesh = {};
    {
        StretchyBuffer<vec3> vertex_positions = {};
        StretchyBuffer<int3> triangle_indices = {};
        {
            FILE *fp = fopen(filename, "r");
            ASSERT(fp);
            char buffer[4096];
            while (fgets(buffer, NELEMS(buffer), fp) != NULL) {
                char prefix[16] = {};
                sscanf(buffer, "%s", prefix);
                if (strcmp(prefix, "f") == 0) {
                    int i, j, k;
                    ASSERT(sscanf(buffer, "%s %d %d %d", prefix, &i, &j, &k) == 4);
                    sbuff_push_back(&triangle_indices, { i - 1, j - 1, k - 1 });
                }
                if (strcmp(prefix, "v") == 0) {
                    real x, y, z;
                    ASSERT(sscanf(buffer, "%s %lf %lf %lf", prefix, &x, &y, &z) == 4);
                    sbuff_push_back(&vertex_positions, { x, y, z });
                }
            }
            fclose(fp);
        }
        // note: don't free the data pointers! (we're stealing them)
        fancy_mesh.num_triangles = triangle_indices.length;
        fancy_mesh.triangle_indices = triangle_indices.data;
        fancy_mesh.num_vertices = vertex_positions.length;
        fancy_mesh.vertex_positions = vertex_positions.data;
    }
    if (transform_vertex_positions_to_double_unit_box) {
        _mesh_utils_transform_vertex_positions_to_double_unit_box(fancy_mesh.num_vertices, fancy_mesh.vertex_positions);
    }
    if (merge_duplicated_vertices) {
        _mesh_utils_indexed_triangle_mesh_merge_duplicated_vertices(&fancy_mesh);
    }
    if (compute_normals) {
        _mesh_utils_indexed_triangle_mesh_alloc_compute_and_store_area_weighted_vertex_normals(&fancy_mesh);
    }
    return fancy_mesh;
}

Soup3D _mesh_utils_soup_TRIANGLES_load(char *filename, bool transform_vertex_positions_to_double_unit_box) {
    Soup3D soup_mesh = {};
    {
        soup_mesh.primitive = SOUP_TRIANGLE_MESH;
        soup_mesh.dimension_of_positions = 3;

        StretchyBuffer<vec3> vertex_positions = {};
        {
            FILE *fp = fopen(filename, "r");
            ASSERT(fp);
            char buffer[4096];
            while (fgets(buffer, NELEMS(buffer), fp) != NULL) {
                real x, y, z;
                ASSERT(sscanf(buffer, "%lf %lf %lf", &x, &y, &z) == 3);
                sbuff_push_back(&vertex_positions, { x, y, z });
            }
            fclose(fp);
        }
        // note: don't free the data pointers! (we're stealing them)
        soup_mesh.num_vertices = vertex_positions.length;
        soup_mesh.vertex_positions = vertex_positions.data;
    }
    if (transform_vertex_positions_to_double_unit_box) {
        _mesh_utils_transform_vertex_positions_to_double_unit_box(soup_mesh.num_vertices, soup_mesh.vertex_positions);
    }
    return soup_mesh;
}
#endif

////////////////////////////////////////////////////////////////////////////////
// #include "utils.cpp"/////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

real utils_random_real(real a = 0, real b = 1) {
    return LERP(real(rand()) / RAND_MAX, a, b);
}

int utils_random_sign() {
    return utils_random_real() < .5 ? -1 : 1;
}

long utils_timestamp_in_milliseconds() { // no promises this is even a little bit accurate
    using namespace std::chrono;
    milliseconds ms = duration_cast<milliseconds>(system_clock::now().time_since_epoch());
    return (long) ms.count();
}


////////////////////////////////////////////////////////////////////////////////
// #include "widget.cpp"////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

bool _widget_drag(real *PV, int num_vertices, real *vertex_positions, real size_in_pixels, real r, real g, real b, real a) {
    if (cow._mouse_owner != COW_MOUSE_OWNER_NONE && cow._mouse_owner != COW_MOUSE_OWNER_WIDGET) return false;
    static real *selected;
    if (selected) { // fornow: allows multiple calls to this function between begin_frame
        bool found = false;
        for (int i = 0; i < num_vertices; ++i) {
            if (selected == vertex_positions + 2 * i) {
                found = true;
                break;
            }
        }
        if (!found) return false;
    }
    real *hot = selected;
    if (!selected) {
        {
            for (int i = 0; i < num_vertices; ++i) {
                real *ptr = vertex_positions + 2 * i;
                real tmp[4] = { ptr[0], ptr[1], 0, 1 };
                _linalg_mat4_times_vec4_persp_divide(tmp, PV, tmp);
                for (int d = 0; d < 2; ++d) tmp[d] -= cow.mouse_position_NDC[d];
                if (_linalg_vecX_squared_length(2, tmp) < pow(.02, 2)) { // todo comparison in pixels
                    hot = ptr;
                }
            }
        }
        if (hot) {
            if (cow.mouse_left_pressed) {
                selected = hot;
            }
        }
    }

    real mouse_xy_World[2];
    _input_get_mouse_position_and_change_in_position_in_world_coordinates(PV, mouse_xy_World, mouse_xy_World + 1, NULL, NULL);
    if (cow.mouse_left_held && selected) {
        for (int d = 0; d < 2; ++d) *(selected + d) = mouse_xy_World[d];
    } else if (cow.mouse_left_released) {
        selected = NULL;
    }

    if (hot || selected) {
        _soup_draw(PV, SOUP_POINTS, _SOUP_XY, _SOUP_RGB, 1, hot, NULL, r, g, b, a, size_in_pixels, false, true);
    }

    cow._mouse_owner = (hot || selected) ? COW_MOUSE_OWNER_WIDGET : COW_MOUSE_OWNER_NONE;
    return selected;
}

#ifdef SNAIL_CPP
template<int D_color = 3> bool widget_drag(mat4 PV, int num_vertices, vec2 *vertex_positions, real size_in_pixels = 0, Vec<D_color> color = { 1.0, 1.0, 1.0 }) {
    STATIC_ASSERT(D_color == 3 || D_color == 4);
    return _widget_drag(PV.data, num_vertices, (real *) vertex_positions, size_in_pixels, color[0], color[1], color[2], D_color == 4 ? color[3] : 1);
}
#endif

////////////////////////////////////////////////////////////////////////////////
// #include "sound.cpp"/////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void _sound_init() {
    ASSERT(cs_init(COW0._window_hwnd__note_this_is_NULL_if_not_on_Windows, 44100, 8192, NULL) == CUTE_SOUND_ERROR_NONE);
    cs_spawn_mix_thread();
    cs_mix_thread_sleep_delay(5);
}

void _sound_reset() {
    cs_stop_all_playing_sounds();
    cs_music_stop(0);
}

void sound_attach_to_gui() {
    real tmp = COW1._sound_music_gui_1_minus_volume;
    gui_slider("music volume                                                                                                                                ", &COW1._sound_music_gui_1_minus_volume, 1.0, 0.0, false);
    if (COW1._sound_music_gui_1_minus_volume != tmp) {
        cs_music_set_volume(1.0f - float(COW1._sound_music_gui_1_minus_volume));
    }
}

int _sound_load(char *filename) {
    ASSERT(COW1._sound_num_loaded < SOUND_MAX_DIFFERENT_FILES);
    ASSERT(strlen(filename) < SOUND_MAX_FILENAME_LENGTH - 1); // ?
    cs_error_t err;
    strcpy(COW1._sound_filenames[COW1._sound_num_loaded], filename);
    COW1._sound_audio_source_ptrs[COW1._sound_num_loaded] = cs_load_wav(filename, &err);
    ASSERT(err == CUTE_SOUND_ERROR_NONE);
    return COW1._sound_num_loaded++;
}

void _sound_play_sound(int audio_source_i) {
    cs_play_sound(COW1._sound_audio_source_ptrs[audio_source_i], cs_sound_params_default());
}

void _sound_loop_music(int audio_source_i) {
    cs_sound_params_t params = cs_sound_params_default();
    params.looped = true;
    cs_music_play(COW1._sound_audio_source_ptrs[audio_source_i], 0);
    cs_music_set_loop(true);
}

int _sound_find_load(char *filename) { // fornow O(n)
    int audio_source_i = -1;
    for (int i = 0; i < COW1._sound_num_loaded; ++i) {
        if (strcmp(COW1._sound_filenames[i], filename) == 0) {
            audio_source_i = i;
            break;
        }
    }
    if (audio_source_i == -1) {
        audio_source_i = _sound_load(filename);
    }
    return audio_source_i;
}

void sound_play_sound(char *filename) {
    _sound_play_sound(_sound_find_load(filename));
}

void sound_loop_music(char *filename) {
    _sound_loop_music(_sound_find_load(filename));
}

////////////////////////////////////////////////////////////////////////////////
// #include "recorder.cpp" /////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#define RECORDER_DAT_FILENAME "codebase/recordings/tmp.dat"
#define MAX_SCRATCH_FILESIZE Gigabytes(25)

void _recorder_draw_pacman(real r, real g, real b, real a, real f) {
    real aspect = window_get_aspect();
    _eso_begin((real *) &cow.I, SOUP_TRIANGLE_FAN, 0.0, false, true);
    real o[2] = { (aspect - .25) / aspect, .75 };
    real radius = .125;
    int N = 32;
    eso_color(r, g, b, a);
    eso_vertex(o[0], o[1]);
    for (int i = 0; i < N; ++i) {
        real theta = f * 2 * PI * i / (N - 1);
        eso_vertex(o[0] + radius * cos(theta) / aspect, o[1] + radius * sin(theta));
    }
    eso_end();
}

void _recorder_begin_frame() { // record
    if ((cow.key_shift_held && cow.key_pressed['`']) || COW0._recorder_out_of_space) {
        COW0._recorder_recording = !COW0._recorder_recording;
        if (COW0._recorder_recording) { // start
            real _width, _height;
            _window_get_size(&_width, &_height);
            COW0._recorder_width = int(_width);
            COW0._recorder_height = int(_height);
            COW0._recorder_size_of_frame = 4 * COW0._recorder_width * COW0._recorder_height;
            COW0._recorder__buffer1 = (byte *) malloc(COW0._recorder_size_of_frame);
            COW0._recorder_buffer2 = (byte *) malloc(COW0._recorder_size_of_frame);
            COW0._recorder_num_frames_recorded = 0;
            COW0._recorder_out_of_space = false;
            char filename[64] = {}; {
                struct tm *timenow;
                time_t now = time(NULL);
                timenow = gmtime(&now);
                strftime(filename, sizeof(filename), "codebase/recordings/%Y-%m-%d--%H-%M-%S.mpg", timenow);
            }
            COW0._recorder_fp_mpg = fopen(filename, "wb");
            if (config.tweaks_record_raw_then_encode_everything_WARNING_USES_A_LOT_OF_DISK_SPACE) {
                COW0._recorder_fp_dat = fopen(RECORDER_DAT_FILENAME, "wb");
            }
        } else { // stop
            if (config.tweaks_record_raw_then_encode_everything_WARNING_USES_A_LOT_OF_DISK_SPACE) {
                ASSERT(COW0._recorder_fp_dat);
                fclose(COW0._recorder_fp_dat);
                COW0._recorder_fp_dat = fopen(RECORDER_DAT_FILENAME, "rb");

                for (int frame = 0; frame < COW0._recorder_num_frames_recorded; ++frame) {
                    fread(COW0._recorder_buffer2, COW0._recorder_size_of_frame, 1, COW0._recorder_fp_dat);
                    jo_write_mpeg(COW0._recorder_fp_mpg, COW0._recorder_buffer2, COW0._recorder_width, COW0._recorder_height, 60);

                    { // progress foo
                        _window_clear_draw_buffer();
                        _recorder_draw_pacman(1.0, 0.0, 0.0, 1.0, 1.0);
                        _recorder_draw_pacman(0.0, 1.0, 0.0, 1.0, LINEAR_REMAP(frame, 0, COW0._recorder_num_frames_recorded - 1, 0.0, 1.0));
                        _window_swap_draw_buffers();
                        glfwPollEvents(); // must poll events to avoid a crash
                    }
                }

                fclose(COW0._recorder_fp_dat);
                ASSERT(remove(RECORDER_DAT_FILENAME) == 0);

                _window_clear_draw_buffer();
            }

            fclose(COW0._recorder_fp_mpg);

            free(COW0._recorder__buffer1);
            free(COW0._recorder_buffer2);
            // COW1.recorder = {};
        }
    }
    if (COW0._recorder_recording) { // save frame
        { // draw mouse
            eso_color(1.0, 1.0, 1.0, 1.0);
            _eso_begin((real *) &cow.NDC_from_Screen, SOUP_TRIANGLE_FAN, 10.0, false, true);
            eso_vertex(cow.mouse_position_Screen[0] + COW0._window_macbook_retina_scale *  0, cow.mouse_position_Screen[1] + COW0._window_macbook_retina_scale *  0);
            eso_vertex(cow.mouse_position_Screen[0] + COW0._window_macbook_retina_scale *  0, cow.mouse_position_Screen[1] + COW0._window_macbook_retina_scale * 14);
            eso_vertex(cow.mouse_position_Screen[0] + COW0._window_macbook_retina_scale *  4, cow.mouse_position_Screen[1] + COW0._window_macbook_retina_scale * 11);
            eso_vertex(cow.mouse_position_Screen[0] + COW0._window_macbook_retina_scale *  6, cow.mouse_position_Screen[1] + COW0._window_macbook_retina_scale * 11);
            eso_vertex(cow.mouse_position_Screen[0] + COW0._window_macbook_retina_scale * 12, cow.mouse_position_Screen[1] + COW0._window_macbook_retina_scale * 11);
            eso_end();
            _eso_begin((real *) &cow.NDC_from_Screen, SOUP_QUADS, 10.0, false, true);
            eso_vertex(cow.mouse_position_Screen[0] + COW0._window_macbook_retina_scale *  4, cow.mouse_position_Screen[1] + COW0._window_macbook_retina_scale * 11);
            eso_vertex(cow.mouse_position_Screen[0] + COW0._window_macbook_retina_scale *  6, cow.mouse_position_Screen[1] + COW0._window_macbook_retina_scale * 11);
            eso_vertex(cow.mouse_position_Screen[0] + COW0._window_macbook_retina_scale *  9, cow.mouse_position_Screen[1] + COW0._window_macbook_retina_scale * 17);
            eso_vertex(cow.mouse_position_Screen[0] + COW0._window_macbook_retina_scale *  7, cow.mouse_position_Screen[1] + COW0._window_macbook_retina_scale * 17);
            eso_end();
        }

        glReadPixels(0, 0, COW0._recorder_width, COW0._recorder_height, GL_RGBA, GL_UNSIGNED_BYTE, COW0._recorder__buffer1);

        if (config.tweaks_record_raw_then_encode_everything_WARNING_USES_A_LOT_OF_DISK_SPACE && ((unsigned long long)(COW0._recorder_num_frames_recorded + 1) * (unsigned long long)(COW0._recorder_size_of_frame)) > MAX_SCRATCH_FILESIZE) {
            COW0._recorder_out_of_space = true;
        } else {
            for (int row = 0; row < COW0._recorder_height; ++row) {
                memcpy(
                        COW0._recorder_buffer2 + row * (4 * COW0._recorder_width),
                        COW0._recorder__buffer1 + (COW0._recorder_height - 1 - row) * (4 * COW0._recorder_width),
                        4 * COW0._recorder_width
                      );
            }

            if (config.tweaks_record_raw_then_encode_everything_WARNING_USES_A_LOT_OF_DISK_SPACE) {
                fwrite(COW0._recorder_buffer2, COW0._recorder_size_of_frame, 1, COW0._recorder_fp_dat);
            } else {
                jo_write_mpeg(COW0._recorder_fp_mpg, COW0._recorder_buffer2, COW0._recorder_width, COW0._recorder_height, 60);
            }

            ++COW0._recorder_num_frames_recorded;
        }

        _recorder_draw_pacman(1.0, 0.0, 0.0, 0.9, 1.0);
    }
}

////////////////////////////////////////////////////////////////////////////////
// #include "meshlib.cpp"///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#include "meshlib.cpp"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// advanced ////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// #include "optimization.cpp"//////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#ifdef SNAIL_CPP
struct SparseMatrixEntry {
    int row;
    int col;
    double val;
};

double *_opt_sparse2dense(int R, int C, int num_entries, SparseMatrixEntry *sparse) {
    double *dense = (double *) calloc(R * C, sizeof(double));
    #define RXC(M, row, col) ((M)[C * (row) + (col)])
    for (int k = 0; k < num_entries; ++k) { RXC(dense, sparse[k].row, sparse[k].col) += sparse[k].val; }
    #undef RXC
    return dense;
}

#ifdef USE_EIGEN
struct EigenTriplet { int row, col; real val; };
void eigenSimplicialCholesky(int N, double *x, int A_length, EigenTriplet *A_data, double *b);
#endif
void opt_solve_sparse_linear_system(int N, double *x, int _A_num_entries, SparseMatrixEntry *_A, double *b) {
    { // checks
        ASSERT(x);
        ASSERT(N);
        ASSERT(_A);
        ASSERT(b);
    }

    if (_A_num_entries == 0) {
        printf("A empty\n");
        memset(x, 0, N * sizeof(double));
        return;
    }

    #ifdef USE_EIGEN
    {
        eigenSimplicialCholesky(N, x, _A_num_entries, (EigenTriplet *) _A, b);
    }
    #else
    {
        do_once { printf("[warn] USE_EIGEN not #define'd; falling back to dense gauss-jordan\n"); };

        // build the augmented matrix
        double *A = _opt_sparse2dense(N, N + 1, _A_num_entries, _A);
        #define NXNP1(M, row, col) ((M)[(N + 1) * (row) + (col)])
        for (int k = 0; k < N; ++k) { NXNP1(A, k, N) = b[k]; }

        { // convert to triangular form (in place)
            // https://en.wikipedia.org/wiki/Gaussian_elimination
            int m = N;
            int n = N + 1;
            int h = 0;
            int k = 0;

            double *scratch = (double *) malloc(n * sizeof(double));
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
                if (ARE_EQUAL(0., NXNP1(A, max_i, k))) {
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
            free(scratch);
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
        free(A);
    }
    #endif
}

void opt_add(double *U, double a) {
    if (U != NULL) {
        *U += a;
    }
};

void opt_add(double *a, int i, vec2 a_i) {
    if (a != NULL) {
        a[2 * i + 0] += a_i[0];
        a[2 * i + 1] += a_i[1];
    }
};

void opt_add(StretchyBuffer<SparseMatrixEntry> *A, int i, int j, mat2 A_ij) {
    if (A != NULL) {
        sbuff_push_back(A, { 2 * i + 0, 2 * j + 0, A_ij(0, 0) });
        sbuff_push_back(A, { 2 * i + 1, 2 * j + 0, A_ij(1, 0) });
        sbuff_push_back(A, { 2 * i + 1, 2 * j + 1, A_ij(1, 1) });
        sbuff_push_back(A, { 2 * i + 0, 2 * j + 1, A_ij(0, 1) });
    }
};

double opt_Vector_dot(int N, double *u, double *v) {
    double ret = 0;
    for (int i = 0; i < N; ++i) {
        ret += u[i] * v[i];
    }
    return ret;
}
#endif

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// top-level functions /////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void _cow_init() {
    ASSERT(!COW0._cow_initialized);

    setvbuf(stdout, NULL, _IONBF, 0); // don't buffer printf

    srand(0);
    // srand((unsigned int) time(NULL));


    _eso_init();
    _window_init();
    _soup_init();
    _itri_init();
    _sound_init();

    COW0._cow_initialized = true;
}

void _cow_reset() {
    ASSERT(COW0._cow_initialized);

    COW1 = {};
    cow._mouse_owner = COW_MOUSE_OWNER_NONE; // fornow

    _sound_reset();
    _window_reset();
}

bool cow_begin_frame() {
    ASSERT(COW0._cow_initialized);

    _window_get_NDC_from_Screen((real *) &cow.NDC_from_Screen);

    { // help force_draw_on_top
        static bool help;
        static bool push_gui_hide_and_disable;
        if (cow.key_shift_held && cow.key_pressed['/']) {
            help = !help;
            if (help) {
                push_gui_hide_and_disable = COW1._gui_hide_and_disable;
            }
            if (!help) {
                COW1._gui_hide_and_disable = push_gui_hide_and_disable;
            }
        }
        // TODO better hiding
        if (help) {
            real box[] = { -1, -1, 1, -1, 1, 1, -1, 1 };
            _soup_draw((real *) &cow.I, SOUP_QUADS, 2, 4, 4, box, NULL, 0.0, 0.0, 0.0, 0.8, 0, false, true);
            COW1._gui_hide_and_disable = false; {
                _gui_begin_frame();
                gui_printf("config.hotkeys_*");
                gui_printf("next app (wraps around) `%s", _gui_hotkey2string(config.hotkeys_app_next));
                gui_printf("previous app (wraps around) `%s", _gui_hotkey2string(config.hotkeys_app_prev));
                gui_printf("quit (next app no wrap) `q");
                gui_printf("quit all `Q");
                gui_printf("exit to main menu `%s", _gui_hotkey2string(config.hotkeys_app_menu));
                gui_printf("display fps counter` \\");
                gui_printf("uncap fps `/");
                gui_printf("(un)hide gui `%s", _gui_hotkey2string(config.hotkeys_gui_hide));
                gui_printf("(un)hide help `?");
                gui_printf("start/stop recording (note: no sound) `~");
                gui_printf("");
                gui_printf("config.tweaks_*");
                gui_checkbox("soup_draw_with_rounded_corners_for_all_line_primitives", &config.tweaks_soup_draw_with_rounded_corners_for_all_line_primitives);
                gui_slider("size_in_pixels_soup_draw_defaults_to_if_you_pass_0_for_size_in_pixels", &config.tweaks_size_in_pixels_soup_draw_defaults_to_if_you_pass_0_for_size_in_pixels, 2.0, 32.0, false);
                gui_checkbox("ASSERT_crashes_the_program_without_you_having_to_press_Enter", &config.tweaks_ASSERT_crashes_the_program_without_you_having_to_press_Enter);
                gui_checkbox("record_raw_then_encode_everything_WARNING_USES_A_LOT_OF_DISK_SPACE", &config.tweaks_record_raw_then_encode_everything_WARNING_USES_A_LOT_OF_DISK_SPACE);
            } COW1._gui_hide_and_disable = true;
        }
    }

    _recorder_begin_frame();
    _window_begin_frame();
    _input_begin_frame();
    _gui_begin_frame();

    { // framerate overlay
        static int measured_fps;
        // request uncapped framerate 
        if (cow.key_pressed['/'] && !cow.key_shift_held) {
            static bool uncapped;
            uncapped = !uncapped;
            glfwSwapInterval(!uncapped);
        }
        { // display fps
            static bool display_fps;
            if (cow.key_pressed['\\']) {
                display_fps = !display_fps;
            }
            if (display_fps) {
                static int fps;
                static long stamp = utils_timestamp_in_milliseconds();
                if (utils_timestamp_in_milliseconds() - stamp > 166) {
                    stamp = utils_timestamp_in_milliseconds();
                    fps = measured_fps;
                    // printf("fps: %d\n", display_fps);
                }
                char text[16] = {};
                snprintf(text, sizeof(text), "fps: %d", fps);
                _text_draw((real *) &cow.NDC_from_Screen, text, 0.0, 0.0, 0.0, (fps < 45) ? 1.0 : 0.0, (fps > 30) ? 1.0 : 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, true);
            }
        }
        real dt;
        // grab and smooth fps
        {
            const int N_MOVING_WINDOW = 5;
            static long prev_stamps[N_MOVING_WINDOW];
            long stamp = utils_timestamp_in_milliseconds();
            measured_fps = (int) round(N_MOVING_WINDOW / (real(stamp - prev_stamps[N_MOVING_WINDOW - 1]) / 1000.));
            for (int i = N_MOVING_WINDOW - 1; i >= 1; --i) {
                prev_stamps[i] = prev_stamps[i - 1];
            }
            prev_stamps[0] = stamp;
            dt = prev_stamps[0] - prev_stamps[1];
        }
        { // cute_sound
            // _sound_update(dt);
        }
    }

    return !(
            cow.key_pressed[config.hotkeys_app_next]
            || cow.key_pressed[config.hotkeys_app_prev]
            || cow.key_pressed['q']
            || cow.key_pressed[config.hotkeys_app_menu]
            );
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// eg ////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#ifdef SNAIL_CPP
void eg_input() {

}

void eg_gl() {

}

void eg_meshlib() {
    Camera3D camera = { 5.0, RAD(0.0) };
    real time = 0.0;
    bool paused = false;
    bool draw_axes = false;

    while (cow_begin_frame()) {
        camera_move(&camera);
        camera_attach_to_gui(&camera);
        gui_checkbox("paused", &paused, 'p');

        mat4 P = camera_get_P(&camera);
        mat4 V = camera_get_V(&camera);
        mat4 PV = P * V;

        mat4 R = RotationAboutY(time);

        mat4 M_wire = Translation(-2.2, 0.0, 0.0) * R;
        mat4 M_smooth = Translation(0.0, 0.0, 0.0) * R;
        mat4 M_matcap = Translation( 2.2, 0.0, 0.0) * R;

        meshlib.soup_bunny.draw(PV * M_wire, monokai.purple);
        meshlib.itri_bunny.draw(P, V, M_smooth, monokai.purple);
        meshlib.itri_bunny.draw(P, cow.I, V * M_matcap, {}, "codebase/matcap.png");

        gui_checkbox("draw_axes", &draw_axes);
        if (draw_axes) {
            meshlib.soup_axes.draw(PV);
        }

        if (!paused) {
            time += .0167;
        }
    }
}

void eg_text() {
    Camera3D camera = { 5.0, RAD(45.0) };
    real time = 0;
    bool paused = false;

    while (cow_begin_frame()) {
        camera_move(&camera);
        mat4 PV = camera_get_PV(&camera);

        text_draw(cow.NDC_from_Screen, "<- mouse", cow.mouse_position_Screen); 

        {
            char *text = "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.";
            int N = int(strlen(text));
            for (int i = 0; i < strlen(text); ++i) {
                char buffer[] = { text[i], '\0' };
                real theta = 5 * LINEAR_REMAP(i, 0, N, 0.0, 2 * PI) - time;
                text_draw(
                        PV,
                        buffer,
                        { cos(theta), LINEAR_REMAP(i, 0, N - 1, 2.5, -2.5), sin(theta) },
                        color_plasma(.5 + .5 * sin(theta))
                        ); 
            }
        }

        if (!paused) {
            time += .0167;
        }
    }
}

void eg_sound() {
    window_set_clear_color(.5 * monokai.blue);

    sound_loop_music("codebase/music.wav");
    while (cow_begin_frame()) {
        if (gui_button("play sound.wav", ' ')) {
            sound_play_sound("codebase/sound.wav");
        }
        sound_attach_to_gui();
    }
}

void eg_soup() {
    Camera2D camera = { 20.0, -4.0 };
    int num_polygon_sides = 16;
    vec2 foo[] = { { -6.0, -6.0 }, { -6.0, 6.0 }, { 6.0, 6.0 }, { 6.0, -6.0 } };
    real size_in_pixels = 12.0;
    bool use_world_units_instead_of_pixels = false;
    bool force_draw_on_top = false;

    while (cow_begin_frame()) {
        camera_move(&camera);
        mat4 PV = camera_get_PV(&camera);

        gui_slider("size_in_pixels", &size_in_pixels, 0, 100, false);
        gui_checkbox("use_world_units_instead_of_pixels", &use_world_units_instead_of_pixels, COW_KEY_TAB);
        gui_checkbox("force_draw_on_top", &force_draw_on_top);
        gui_slider("num_polygon_sides", &num_polygon_sides, 0, 32, 'j', 'k', false);

        eso_begin(PV, SOUP_LINE_LOOP, size_in_pixels, use_world_units_instead_of_pixels, force_draw_on_top); {
            eso_color(monokai.green);
            for (int i = 0; i < num_polygon_sides; ++i) {
                real theta = real(i) / real(num_polygon_sides) * 2.0 * PI;
                real r = 6.0 * sqrt(2.0);
                eso_vertex(r * V2(cos(theta), sin(theta)));
            }
        } eso_end();

        widget_drag(PV, 4, foo, size_in_pixels, monokai.yellow);
        soup_draw(PV, SOUP_QUADS, 4, foo, NULL, V4(monokai.red, .5), 0, use_world_units_instead_of_pixels, force_draw_on_top);

        vec2 s_mouse = input_get_mouse_position_in_world_coordinates(PV);
        eso_begin(PV, SOUP_POINTS, size_in_pixels, use_world_units_instead_of_pixels, force_draw_on_top);
        eso_color(monokai.blue);
        eso_vertex(s_mouse);
        eso_end();
    }
}

void eg_kitchen_sink() {
    Camera3D camera = { 5.0, RAD(0.0) };
    real time = 0.0;
    bool paused = false;
    bool draw_axes = false;

    StretchyBuffer<vec2> trace = {};

    sound_loop_music("codebase/music.wav");
    while (cow_begin_frame()) {
        camera_move(&camera);
        camera_attach_to_gui(&camera);

        if (gui_button("play sound.wav", ' ')) {
            sound_play_sound("codebase/sound.wav");
        }
        sound_attach_to_gui();

        // todo paint trails coming off of mouse

        gui_checkbox("paused", &paused, 'p');

        mat4 P = camera_get_P(&camera);
        mat4 V = camera_get_V(&camera);
        mat4 PV = P * V;

        mat4 R = RotationAboutY(time);

        mat4 M_wire = Translation(-2.2, 0.0, 0.0) * R;
        mat4 M_smooth = Translation(0.0, 0.0, 0.0) * R;
        mat4 M_matcap = Translation( 2.2, 0.0, 0.0) * R;

        vec3 color = !(cow.mouse_left_held && !cow._mouse_owner) ? monokai.purple : monokai.blue;
        meshlib.soup_bunny.draw(PV * M_wire, color);
        meshlib.itri_bunny.draw(P, V, M_smooth, color);
        meshlib.itri_bunny.draw(P, cow.I, V * M_matcap, {}, "codebase/matcap.png");

        gui_checkbox("draw_axes", &draw_axes, COW_KEY_TAB);
        if (draw_axes) {
            meshlib.soup_axes.draw(PV);
        }

        vec2 s_mouse = cow.mouse_position_NDC;
        if (trace.length == 0 || squaredNorm(trace[trace.length - 1] - s_mouse) > .0001) {
            sbuff_push_back(&trace, s_mouse);
        }
        soup_draw(cow.I, SOUP_LINE_STRIP, trace.length, trace.data, NULL, color_rainbow_swirl(time), 0, false, true);
        if (gui_button("clear trace", 'r')) {
            sbuff_free(&trace);
        }

        text_draw(cow.NDC_from_Screen, "  :3", cow.mouse_position_Screen); 

        {
            char *text = "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.";
            int N = int(strlen(text));
            for (int i = 0; i < strlen(text); ++i) {
                char buffer[] = { text[i], '\0' };
                real theta = time - 5.0 * LINEAR_REMAP(i, 0, N, 0.0, 2 * PI);
                text_draw(
                        PV,
                        buffer,
                        { cos(theta), LINEAR_REMAP(i, 0, N - 1, 2.5, -2.5), sin(theta) },
                        color_plasma(.5 + .5 * sin(theta))
                        ); 
            }
        }

        if (!paused) {
            time += .0167;
        }

        {
            vec3 o, x, y, z;

            o = { 0.0, 2.0, 0.0 };
            x = { 0.1, 0.0, 0.0 };
            y = { 0.0, 0.1, 0.0 };
            z = { 0.0, 0.0, 0.1 };

            eso_begin(PV, SOUP_LINES);
            eso_color(monokai.red);
            eso_vertex(o);
            eso_vertex(o + x);
            eso_color(monokai.green);
            eso_vertex(o);
            eso_vertex(o + y);
            eso_color(monokai.blue);
            eso_vertex(o);
            eso_vertex(o + z);
            eso_end();
        }
    }
}
#endif

void _eg_no_snail() {
    Camera2D camera = { 5.0, -1.0 };
    int num_polygon_sides = 16;
    real foo[] = { -1.5, -1.5, -1.5, 1.5, 1.5, 1.5, 1.5, -1.5 };
    real size_in_pixels = 3.0;
    real PV[16];

    while (cow_begin_frame()) {
        camera_move(&camera);
        _camera_get_PV(&camera, PV);

        gui_slider("num_polygon_sides", &num_polygon_sides, 0, 32, 'j', 'k', false);
        gui_slider("size_in_pixels", &size_in_pixels, 0.0, 10.0, false);

        _eso_begin(PV, SOUP_LINE_LOOP, size_in_pixels, false, false); {
            eso_color(0.0, 1.0, 0.0);
            for (int i = 0; i < num_polygon_sides; ++i) {
                real theta = real(i) / real(num_polygon_sides) * 2.0 * PI;
                real r = 1.5 * sqrt(2.0);
                eso_vertex(r * cos(theta), r * sin(theta));
            }
        } eso_end();

        _widget_drag(PV, 4, foo, size_in_pixels, 1.0, 1.0, 0.0, 1.0);
        _soup_draw(PV, SOUP_QUADS, _SOUP_XY, _SOUP_RGB, 4, foo, NULL, 1.0, 0.0, 0.0, 1.0, 0, false, false);

        real s_mouse[2];
        _input_get_mouse_position_and_change_in_position_in_world_coordinates(
                PV,
                s_mouse,
                s_mouse + 1,
                NULL,
                NULL
                );
        _eso_begin(PV, SOUP_POINTS, size_in_pixels, false, false);
        eso_color(0.0, 0.0, 1.0);
        eso_vertex(s_mouse[0], s_mouse[1]);
        eso_end();
    }
}

#ifdef SNAIL_CPP
#define _APP_EXAMPLES_ALL() \
    APP(_eg_no_snail); \
    APP(eg_soup); \
    APP(eg_meshlib); \
    APP(eg_text); \
    APP(eg_sound);
#else 
#define _APP_EXAMPLES_ALL() APP(_eg_no_snail);
#endif

#endif
