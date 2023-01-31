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

int main() {
    APPS {
        APP(eg_kitchen_sink);
        APP(app_hello_snail);
        APP(app_hello_cow);
        APP(app_hello_extra_credit);
    }

    return 0;
}
