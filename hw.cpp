// !! This is the hw.cpp for hw1.                                               
//                                                                              
// If you cloned today expecting to start working on hw0, please...             
//                                                                              
// - open the file explorer / finder                                            
// - navigate into folder old_hw\hw0                                            
// - copy hw.cpp                                                                
// - navigate back to the root folder (CSCI-371)                                
// - paste hw.cpp inside, replacing the file that's already there (this file)   
//                                                                              
// You should then be good to go :)                                             
//                                                                              
// Please ask me ASAP if you have any questions.                                


// if you find a bug in the hw, please report it on Github Issues or email      
//                                        the class will be rewarded with snacks

#define _CRT_SECURE_NO_WARNINGS
// #define COW_CRASH_ON_FLOATING_POINT_EXCEPTIONS
#include "snail.cpp"
#include "cow.cpp"
#include <vector>
#include <iostream>


//                                                                              
// documentation                                                                
//                                                                              

#if 0
// usage code
#endif

// you can quit the current app by pressing the Q key                           

// press the \ key to display an fps counter                                    
// (macbook users see README.txt if you don't have a 60 fps cap)                
// press the / key to uncap the framerate                                       

// here is the actual basic_draw function signature (at least for this week)    
// if vertex_colors is NULL, then it draws all vertices with the same color,    
// specifically { r_fallback, g_fallback, b_fallback, a_fallback }              
// if overlay is true, then we draw on top of _everything_                      
//                                                                              
// void basic_draw(                                                             
//         int primitive,                                                       
//         double *transform,                                                   
//         int dimension_of_positions,                                          
//         int dimension_of_colors,                                             
//         int num_vertices,                                                    
//         double *vertex_positions,                                            
//         double *vertex_colors = NULL,                                        
//         double r_fallback = 1,                                               
//         double g_fallback = 1,                                               
//         double b_fallback = 1,                                               
//         double a_fallback = 1,                                               
//         double size_in_pixels = 0,                                           
//         bool overlay = false,                                                
//         double r_wireframe = 1,                                              
//         double g_wireframe = 1,                                              
//         double b_wireframe = 1,                                              
//         double a_wireframe = 1);                                             
//                                                                              
// snail gives us a couple convenient wrappers, which we can think of as        
//                                                                              
//   // single color                                                            
//   void basic_draw(                                                           
//           int primtive,                                                      
//           mat4 transform,                                                    
//           int num_vertices,                                                  
//           vecX *vertex_positions,                                            
//           vec3 color,                                                        
//           double size_in_pixels = 0,                                         
//           bool overlay = false);                                             
//                                                                              
//   // per-vertex color                                                        
//   void basic_draw(                                                           
//           int primtive,                                                      
//           mat4 transform,                                                    
//           int num_vertices,                                                  
//           vecX *vertex_positions,                                            
//           vec3 *vertex_colors,                                               
//           double size_in_pixels = 0,                                         
//           bool overlay = false);                                             
//                                                                              
// here vecX means you can use vec2 or vec3 (cow.cpp uses templates for this)   

#if 0
static Camera2D camera = { 5 };
mat4 PV = camera_get_PV(&camera);
vec2 foo[3] = { { 0, 0 }, { 1, 0}, { 0, 1} };
basic_draw(TRIANGLES, PV, 3, foo, monokai.white);
#endif

// this function lets you drag 2D vertices                                      
// the arguments with default arguments specify the size and color of the dot   
// that pops up when you hover over a point                                     
//                                                                              
// bool widget_drag(                                                            
//         double *PV,                                                          
//         int num_vertices,                                                    
//         double *vertex_positions,                                            
//         double size_in_pixels = 0,                                           
//         double r = 1,                                                        
//         double g = 1,                                                        
//         double b = 1,                                                        
//         double a = 1);                                                       
//                                                                              
// if we include snail we can call this wrapper if we prefer                    
//                                                                              
// bool widget_drag(                                                            
//         mat4 PV,                                                             
//         int num_vertices,                                                    
//         vec2 *vertex_positions,                                              
//         double size_in_pixels = 0,                                           
//         vec3 color = monokai.white);                                         

#if 0
static Camera2D camera = { 5 };
mat4 PV = camera_get_PV(&camera);
static vec2 foo[3] = { { 0, 0 }, { 1, 0 }, { 0, 1 } };
basic_draw(LINE_LOOP, PV, 3, foo, monokai.white);
widget_drag(PV, 3, foo);
#endif

// imgui_slider lets us scrub an int t between bounds a and b.                  
// pressing the key j will decrement t, pressing the key k will increment it    
// if loop is true, e.g. going above b will take us back to a                   
//                                                                              
// void imgui_slider(                                                           
//         char *name,                                                          
//         int *t,                                                              
//         int a,                                                               
//         int b,                                                               
//         char j = 0,                                                          
//         char k = 0,                                                          
//         bool loop = false);                                                  
//                                                                              
// there is also a version for doubles                                          
//                                                                              
// void imgui_slider(char *name, double *t, double a, double b);                

#if 0
static int foo = 0;
imgui_slider("foo", &foo, 0, 100, 'j', 'k', true);
#endif

// here is how we can access user input                                         
//                                                                              
// input.key_pressed[...]                                                       
// input.key_released[...]                                                      
// input.key_held[...]                                                          
// input.key_toggle[...]                                                        
// input.mouse_left_pressed                                                     
// input.mouse_left_held                                                        
// input.mouse_left_released                                                    
// input.mouse_right_pressed                                                    
// input.mouse_right_held                                                       
// input.mouse_right_released                                                   
//                                                                              
// void input_get_mouse_position_and_change_in_position_in_world_coordinates(   
//         double *PV,                                                          
//         double *mouse_x_world,                                               
//         double *mouse_y_world,                                               
//         double *mouse_dx_world = NULL,                                       
//         double *mouse_dy_world = NULL);                                      
//                                                                              
// which has these snail wrappers                                               
//                                                                              
// vec2 input_get_mouse_position_in_world_coordinates(mat4 PV);                 
// vec2 input_get_mouse_change_in_position_in_world_coordinates(mat4 PV);       

#if 0
if (input.key_pressed['a']) { printf("you pressed the A key!\n"); }

static Camera2D camera = { 5 };
camera_move(&camera);
mat4 PV = camera_get_PV(&camera);
if (input.mouse_left_pressed) {
    vec2 s_mouse = input_get_mouse_position_in_world_coordinates(PV);
    printf("you lefted clicked at (%lf, %lf) in world coordinates\n", s_mouse.x, s_mouse.y);
}
#endif

// NELEMS(fixed_size_array) gives the number of elements in a fixed-size array  
// use at your own risk; careful not to call it on a pointer                    


//                                                                              
// implementation                                                               
//                                                                              

// begin submission                                                             

void hw1a() {
    while (begin_frame()) {
        static Camera2D camera = { 5 };
        camera_move(&camera);
        mat4 PV = camera_get_PV(&camera);

        // todo draggable red square outline

        // todo green circle outline with slider to change number of vertices
    }
}

void hw1b() {
    // Camera3D camera = { 5, RAD(45) };
}

void hw1c() {
    #if 0
    StretchyBuffer buffer = {};

    ASSERT(buffer.length == 0);
    ASSERT(buffer.capacity == 0);

    sbuff_push_back(&buffer, {});

    ASSERT(buffer.length == 1);
    ASSERT(buffer.capacity == 16);

    int N = 2048;
    for (int i = 0; i < N; ++i) {
        double theta = double(i) / 64 * 2 * PI;
        double r = double(i) / N;
        sbuff_push_back(&buffer, r * V2(cos(theta), sin(theta)));
    }

    while (begin_frame()) {
        static Camera2D camera = { .8 };
        camera_move(&camera);
        mat4 PV = camera_get_PV(&camera);

        static double time = 0;
        time += .0167;

        basic_draw(LINE_STRIP, PV * RotationZ(-10 * time), buffer.length, buffer.data, color_rainbow_swirl(.1 * time));
    }

    // note: we don't technically have to free here if the program is going
    // to exit anyway; i just needed to call free somewhere to test it :)  
    sbuff_free(&buffer);
    ASSERT(buffer.length == 0);
    ASSERT(buffer.capacity == 0);
    ASSERT(buffer.data == NULL);
    #endif
}


void hw1d() {
}



void hw() {
    hw1a();
    hw1b();
    hw1c();
    hw1d();

    // you can add additional apps to demo extra credit here
}

// end submission                                                               


int main() {
    init(false, "hw1 :)");
    hw();
    return 0;
}


