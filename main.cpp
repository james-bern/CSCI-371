// #define COW_NO_SOUND
#define JIM_IS_JIM
#include "include.cpp"

#include "wrl/space_fish.cpp"
#include "wrl/interval_timer.cpp"




int main() {
    {
        {
            config.tweaks_record_raw_then_encode_everything_WARNING_USES_A_LOT_OF_DISK_SPACE = true;
            config.tweaks_soup_draw_with_rounded_corners_for_all_line_primitives = true;
        }

        APPS {
            // APP(eg_soup);
            // APP(eg_texture);
            // APP(app_diego);
            // APP(eg_kitchen_sink);
            // APP(app_interval_timer);
            APP(app_space_fish_2D);
            APP(_eg_no_snail);
            // APP(app_space_fish_exploration);
            // APP(app_space_fish_arduino);
        }
    }
    return 0;
}
