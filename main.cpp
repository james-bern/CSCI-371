// TODO sol.cpp solution to all the homeworks from last semester (will tell you exactly what functions to document)
// doughnut 

// TODO user just calls functions (no need to interact with state)
// TODO make it one file
// TODO make it FLAT (prefixes for everything) -- no anonymous structs

#include "include.cpp"


// #include "10sol.cpp"

// new features
// ------------
// ? TRIANGLE_MESH and QUAD_MESH done in one pass with custom shaders
// text entry box for int, real, char * (then have bool to optionally add to slider)
// drop down menu for enums?
// fix transparent window stuff (blend mode?--background shouldn't show through)

// old features to modernize / port
// ------------
// fancy_draw  

// starter code
// ------------
// custom_shader

// homework
// ------------
// matcap
// raymarching -- vec4 (alpha)




int main() {
    config.tweaks_record_raw_then_encode_everything_WARNING_USES_A_LOT_OF_DISK_SPACE = true;
    config.tweaks_soup_draw_with_rounded_corners_for_all_line_primitives = true;

    APPS {
        // APP(exam10);
        APP(eg_soup);
        APP(eg_kitchen_sink);
        // _APP_EXAMPLES_ALL();
    }

    return 0;
}
