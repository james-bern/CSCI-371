#include "include.cpp"

void app_treasure() {
    StretchyBuffer<int> numbers = {};
    sbuff_push_back(&numbers, 4);
    sbuff_push_back(&numbers, 2);
    ASSERT(numbers.length == 2);
    printf("%d", numbers[0]);
    printf("%d", numbers[1]);
    sbuff_free(&numbers);

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
        // TODO add this line to wiki
        soup_draw(PV, SOUP_POINTS, vertex_positions.length, vertex_positions.data, NULL, monokai.red);

    }
}


int main() {
    config.tweaks_record_raw_then_encode_everything_WARNING_USES_A_LOT_OF_DISK_SPACE = true;
    config.tweaks_soup_draw_with_rounded_corners_for_all_line_primitives = true;

    APPS {
        APP(app_treasure);
        APP(eg_kitchen_sink);
    }

    return 0;
}
