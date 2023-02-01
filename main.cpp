#include "include.cpp"

void hw00a() {
    Camera2D camera = { 3.0 };

    vec3 color = V3(1.0, 0.0, 1.0);
    vec2 vertex_positions[3] = {
        V2(0.0, 0.0),
        V2(1.0, 0.0),
        V2(0.0, 1.0),
    };

    vec2 edge_1 = vertex_positions[1] - vertex_positions[0];
    vec2 edge_2 = vertex_positions[2] - vertex_positions[0];
    real triangle_area = .5 * ABS(cross(edge_1, edge_2));
    printf("the triangle's area is %lf\n", triangle_area);



    while (cow_begin_frame()) {
        camera_move(&camera);
        mat4 PV = camera_get_PV(&camera);

        soup_draw(PV, SOUP_TRIANGLES, 3, vertex_positions, NULL, color);
    }
}

int main() {
    APPS {
        APP(hw00a);
    }
    return 0;
}


