#include "include.cpp"

void hw00() {
    while (cow_begin_frame()) {

    }
}

int main() {
    APPS {
        APP(eg_kitchen_sink);
        APP(hw00);
    }
    return 0;
}

#if 0
// #define COW_NO_STYLE_GUIDE
#include "include.cpp"

void hw00a() {
    vec2 s = {};

    Camera2D camera = { 3.0 };
    while(cow_begin_frame()) {
        camera_move(&camera);
        mat4 PV = camera_get_PV(&camera);

        if (gui_button("randomize!")) {
            s = { random_real(-1.0, 1.0), random_real(-1.0, 1.0) };
        }

        eso_begin(PV, SOUP_LINE_LOOP); {
            eso_color(monokai.white);
            eso_vertex(-1.0, -1.0);
            eso_vertex(-1.0,  1.0);
            eso_vertex( 1.0,  1.0);
            eso_vertex( 1.0, -1.0);
        } eso_end();

        widget_drag(PV, 1, &s);
        soup_draw(PV, SOUP_POINTS, 1, &s, NULL, monokai.purple);
    }
}

int main() {
    APPS {
        APP(eg_kitchen_sink);
        APP(hw00a);
    }

    return 0;
}
#endif
