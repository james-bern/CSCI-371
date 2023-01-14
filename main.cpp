#include "include.cpp"

void app_treasure() {
    StretchyBuffer<vec3> vertex_positions = {}; {
        FILE *fp = fopen("TreasureChest.obj", "r");
        ASSERT(fp);
        char line[4096];
        while (fgets(line, _COUNT_OF(line), fp) != NULL) {
            char prefix[64] = {};
            sscanf(line, "%s", prefix);
            if (strcmp(prefix, "v") == 0) {
                real x, y, z;
                ASSERT(sscanf(line, "%s %lf %lf %lf", prefix, &x, &y, &z) == 4);
                sbuff_push_back(&vertex_positions, { x, y, z });
            }
        }
        fclose(fp);
    }

    Camera3D camera = { 10.0, RAD(0) };
    while (cow_begin_frame()) {
        camera_move(&camera);
        mat4 PV = camera_get_PV(&camera);
        vec2 s_mouse = mouse_get_position(PV);
        soup_draw(PV, SOUP_POINTS, 1, &s_mouse, NULL, monokai.red);
        // TODO add this line to wiki
        soup_draw(PV, SOUP_POINTS, vertex_positions.length, vertex_positions.data, NULL, monokai.red);
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
        APP(app_treasure);
        APP(eg_kitchen_sink);
    }

    return 0;
}
