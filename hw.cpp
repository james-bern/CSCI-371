#define COW_CRASH_ON_FLOATING_POINT_EXCEPTIONS

#define _CRT_SECURE_NO_WARNINGS
#include "snail.cpp"
#include "cow.cpp"

// // documentation ////////////////////////////////////////////////////////////
//                                                                              
// // random double in interval [a, b]                                          
// double util_random_double(double a = 0, double b = 1);                       
//                                                                              
// double norm(vecX v); // length of a vector v                                 
//                                                                              
// matX inverse(matX M); // inverse of matrix M                                 
//                                                                              
// // useful macros (see hw2 writeup for explanation)                           
// #define LERP(t, a, b)          ((1 - (t)) * (a) + (t) * (b))                 
// #define INVERSE_LERP(c, a, b)  (((c) - (a)) / double((b) - (a)))             
// #define CLAMP(t, a, b)         MIN(MAX(t, a), b)                             
// #define CLAMPED_LERP(t, a, b)  LERP(CLAMP(t, 0, 1), a, b)                    
// #define COS_LERP(t, a, b)      LERP(.5 - .5 * cos((t)*PI), a, b)             


// begin submission                                                             

void hw2a() {
    init();
    // todo: calculate and print p
    // todo: calculate and print distance from p to q
}

void hw2b() {
    init();
    Camera2D camera = { 10 };
    vec2 circle_center = V2(0, 0);
    double circle_radius = 2;
    vec2 test_point = V2(0, 0);
    while (begin_frame()) {
        camera_move(&camera);
        mat4 PV = camera_get_PV(&camera);
        imgui_slider("r", &circle_radius, 0, 5);
        imgui_slider("x", &circle_center.x, -5, 5);
        imgui_slider("y", &circle_center.y, -5, 5);

        bool inside = false; // todo: set to the correct value

        basic_draw(POINTS, PV, 1, &test_point, (inside) ? monokai.green : monokai.red);
        widget_drag(PV, 1, &test_point, 0, V4(1, 1, 1, .5));

        gl_PV(PV);
        gl_begin(LINE_LOOP);
        gl_color(monokai.yellow);
        for (double theta = 0; theta < 2 * PI; theta += .1) {
            // gl_vertex(circle_center.x + circle_radius * cos(theta), circle_center.y + circle_radius * sin(theta));
            // gl_vertex(circle_center + circle_radius * V2(cos(theta), sin(theta)));
            gl_vertex(circle_center + circle_radius * e_theta(theta));
        }
        gl_end();
    }
}

void hw2c() {
    init();
    // todo
}

void hw2d() {
    init();
    int the_first_frame_that_k_ends_in_47 = 0; // todo: set to the correct value

    // begin don't modify this code
    if (!the_first_frame_that_k_ends_in_47) {
        unsigned int k = 0;
        int frame = 0;
        while (begin_frame()) {
            unsigned int m_w = 1 + k;
            unsigned int m_z = 2 + k;
            m_z = 36969 * (m_z & 65535) + (m_z >> 16);
            m_w = 18000 * (m_w & 65535) + (m_w >> 16);
            k += (m_z << 17) + m_w;
            ++frame;
        }
    } else {
        while (begin_frame()) {
            char r = char(the_first_frame_that_k_ends_in_47 * 233.111111111);
            char g = char(the_first_frame_that_k_ends_in_47 * 2833.33333333);
            g = char(g + 1);
            char b = char(the_first_frame_that_k_ends_in_47 * 210.666666667);
            clear_draw_buffer(r / 255., g / 255., b / 255., 1);
        }
    }
    // end don't modify this code
}

void hw2e() {
    init();
    Camera2D camera = { 5 };
    bool playing = false;
    double t = 0;
    double a = 0;
    double b = 2;
    while (begin_frame()) {
        camera_move(&camera);
        mat4 PV = camera_get_PV(&camera);

        imgui_checkbox("playing", &playing, 'p');
        if (playing) {
            t += 1. / 60;
        }
        if (imgui_button("reset", 'r')) {
            t = 0;
        }
        imgui_slider("t", &t, -2, 2);

        imgui_slider("a", &a, -4, 4);
        imgui_slider("b", &b, -4, 4);

        gl_PV(PV);
        gl_begin(LINES);
        gl_color(monokai.blue);
        gl_vertex(a, 10);
        gl_vertex(a, -10);
        gl_color(monokai.orange);
        gl_vertex(b, 10);
        gl_vertex(b, -10);
        gl_end();

        // todo: draw dots
    }
}

void hw2f() {
    init();
    Camera2D camera = { 3 };
    vec2 a = V2(-1, -1);
    vec2 b = V2(1, -1);
    vec2 c = V2(0, .5 * sqrt(3));
    vec2 p = V2(0, 0);
    vec3 alpha_beta_gamma = V3(1, 1, 1); // (alpha, beta, gamma)
    while (begin_frame()) {
        camera_move(&camera);
        mat4 PV = camera_get_PV(&camera);

        // todo: calculate alpha_beta_gamma

        vec3 color = alpha_beta_gamma.x * V3(1, 0, 0) + alpha_beta_gamma.y * V3(0, 1, 0) + alpha_beta_gamma.z * V3(0, 0, 1);
        basic_draw(POINTS, PV, 1, &p, color);
        widget_drag(PV, 1, &p, 15, color);

        imgui_readout("alpha", &alpha_beta_gamma.x);
        imgui_readout("beta", &alpha_beta_gamma.y);
        imgui_readout("gamma", &alpha_beta_gamma.z);

        gl_PV(PV);
        gl_begin(LINE_LOOP);
        gl_color(V3(1, 0, 0)); gl_vertex(a);
        gl_color(V3(0, 1, 0)); gl_vertex(b);
        gl_color(V3(0, 0, 1)); gl_vertex(c);
        gl_end();
        widget_drag(PV, 1, &a);
        widget_drag(PV, 1, &b);
        widget_drag(PV, 1, &c);

    }
}

void hw2g() {
    // todo
}

void hw() {
    hw2a();
    hw2b();
    hw2c();
    hw2d();
    hw2e();
    hw2f();
    hw2g();
}

// end submission                                                               

int main() {
    // please change this call to init(...)        
    // so you like where the popup window spawns :)
    //                                             
    // void init(                                  
    //     bool transparent_framebuffer,           
    //     char *window_title,                     
    //     int screen_height_in_pixels,            
    //     int window_top_left_init_x_in_pixels,   
    //     int window_top_left_init_y_in_pixels    
    // );                                          
    init(false, "hw2", 540, 0, 100);
    hw();
    return 0;
}
