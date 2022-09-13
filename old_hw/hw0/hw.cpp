// hw guidelines                                                                
// - follow all policies listed here: https://james-bern.github.io/csci-371     
// - don't delete the begin and end submission comments                         
// - put your code in between the begin and end submission comments             
// - otherwise basically anything goes; you will be graded on output, not style 
//   - name things whatever you like; use any number of spaces and tabs;        
//     define new functions, structs, etc.; use global an static variables,     
//     lambda functions, macros*, etc.;                                         
//               *please #undef your #define's before the end submission comment

#define _CRT_SECURE_NO_WARNINGS
#include "snail.cpp"
#include "cow.cpp"
// ^ feel free to experiment with these files; copy parts into hw.cpp; etc.     
// however, your submitted hw.cpp must work with the cow and snail.cpp i shipped
// PS if you find a bug, you can report it on Github Issues*                    
//                                            *you will be rewarded with snacks 
//                                                                              
// v feel free to include any C/C++11 standard library, or anything in ext/stb  
#include <vector>
#include <iostream>


// documentation                                                                
// ---------------------------------------------------------------------------- 
// gl_begin(int primitive double size_in_pixels);       starts drawing          
// gl_color(double r, double g, double b double a = 1); sets the current color  
// gl_vertex(double x, double y, double z = 0);         pushes a vertex         
// gl_end();                                            draws the primitive     
//                                                                              
// - useful primitives: POINTS, LINES, LINE_STRIP, LINE_LOOP, TRIANGLES, QUADS  


// begin submission                                                             

void hw() {
    while (begin_frame()) {
        static Camera2D camera = { 5 };
        camera_move(&camera);
        mat4 PV = camera_get_PV(&camera);

        gl_PV(PV);
        gl_begin(LINE_LOOP);
        gl_color(1, 1, 1);
        gl_vertex(0, 0);
        gl_vertex(1, 0);
        gl_vertex(0, 1);
        gl_end();
    }                                                              
}

// end submission                                                               


int main() {
    init();
    hw();
    return 0;
}


