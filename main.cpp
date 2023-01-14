#include "include.cpp"

void app_sketch() {
    Camera3D camera = { 3.0, RAD(0) };
    while (cow_begin_frame()) {
        camera_move(&camera);
        mat4 PV = camera_get_PV(&camera);
        meshlib.soup_teapot.draw(PV, monokai.green);
    }
}

int main() {
    config.tweaks_record_raw_then_encode_everything_WARNING_USES_A_LOT_OF_DISK_SPACE = true;
    config.tweaks_soup_draw_with_rounded_corners_for_all_line_primitives = true;

    vec2 a = { 1.0, 2.0 };
    vec2 b = a + V2(5.0, 3.0);

    APPS {
        APP(eg_shader);
        // APP(eg_texture);
        // APP(app_sketch);
        APP(eg_kitchen_sink);
    }

    return 0;
}
