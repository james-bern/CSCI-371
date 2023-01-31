#include "include.cpp"

void app_hello_snail() {
    while(cow_begin_frame()) {

    }
}

void app_hello_cow() {
    while(cow_begin_frame()) {

    }
}

void app_hello_extra_credit() {
    while(cow_begin_frame()) {

    }
}

// todo code flow setup/loop
// todo snail

// todo random_vec2(double a, double b)
// todo keyboard intpu

int main() {
    static int i = 0;
    APPS {
        APP(eg_kitchen_sink);
        APP(app_hello_snail);
        APP(app_hello_cow);
        APP(app_hello_extra_credit);
    }

    return 0;
}
