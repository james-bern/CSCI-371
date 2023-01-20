#define JIM_IS_JIM
#include "include.cpp"

int main() {
    {
        config.tweaks_record_raw_then_encode_everything_WARNING_USES_A_LOT_OF_DISK_SPACE = true;
        config.tweaks_soup_draw_with_rounded_corners_for_all_line_primitives = true;
    }

    APPS {
        APP(eg_soup);
        APP(eg_texture);
        APP(eg_kitchen_sink);
    }

    return 0;
}
