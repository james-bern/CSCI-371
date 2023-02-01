// #define COW_NO_STYLE_GUIDE
#include "include.cpp"

void hw00() {
    while(cow_begin_frame()) {

    }
}

int main() {
    APPS {
        APP(eg_kitchen_sink);
        APP(hw00);
    }

    return 0;
}
