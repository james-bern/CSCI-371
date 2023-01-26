#define JIM_IS_JIM
#include "include.cpp"
#include "diego.cpp"

// a b c d e 0 0 0 0 0 0 0 0 
// 0 1 2 3 4 5

template <typename T> struct Vector {
    int length; // how_many_elements_are_currently_stored
    int capacity; // how_many_elements_can_we_currently_fit
    T *data;

    template <typename T> void push_back(T val) {
        if (capacity == length) {
            capacity = (capacity == 0) ? 16 : 2 * capacity;
            data = (T *) realloc(data, capacity * sizeof(T));
        }
        data[length++] = val;
    }
};

int main() {
    {
        // {
        //     config.tweaks_record_raw_then_encode_everything_WARNING_USES_A_LOT_OF_DISK_SPACE = true;
        //     config.tweaks_soup_draw_with_rounded_corners_for_all_line_primitives = true;
        // }

        // APPS {
        //     // APP(eg_soup);
        //     // APP(eg_texture);
        //     APP(app_diego);
        //     APP(eg_kitchen_sink);
        // }
    }

    Vector<double> foo = Vector<double>();
    foo.push_back(1.0);
    foo.push_back(42.0);
    foo.push_back(3.14);
    double bar = 1.0;
    while (1) {
        foo.push_back(bar);
        bar += 1.0;
    }

    return 0;
}
