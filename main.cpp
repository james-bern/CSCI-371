#include "include.cpp"

real get_triangle_signed_area(vec2 *vertex_positions) {
    // vertex_positions should be a pointer to three contiguous vec2's
    // convention: clockwise is positive
    vec2 edge_1 = vertex_positions[1] - vertex_positions[0];
    vec2 edge_2 = vertex_positions[2] - vertex_positions[0];
    return .5 * cross(edge_1, edge_2);
}

void hw00a() {
    Camera2D camera = { 3.0 };

    vec3 color = V3(1.0, 0.0, 1.0);
    vec2 vertex_positions[3] = {
        V2(0.0, 0.0),
        V2(1.0, 0.0),
        V2(0.0, 1.0),
    };

    {
        real signed_area = get_triangle_signed_area(vertex_positions);
        printf("area_0 %lf\n", signed_area);
   }

    while (cow_begin_frame()) {
        camera_move(&camera);
        mat4 PV = camera_get_PV(&camera);

        real signed_area = get_triangle_signed_area(vertex_positions);
        gui_readout("signed_area", &signed_area);

        widget_drag(PV, 3, vertex_positions);
        soup_draw(PV, SOUP_TRIANGLES, 3, vertex_positions, NULL, color);
    }
}


int main() {
    APPS {
        APP(hw00a);
    }
    return 0;
}


