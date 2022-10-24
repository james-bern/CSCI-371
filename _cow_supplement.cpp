// stuff you already implemented in the homeworks :)

template <typename T> struct StretchyBuffer {
    int length;
    int capacity;
    T *data;
    T &operator [](int index) { return data[index]; }
};

template <typename T> void sbuff_push_back(StretchyBuffer<T> *buffer, T element) {
    if (buffer->capacity == 0) {
        buffer->capacity = 16;
        buffer->data = (T *) malloc(buffer->capacity * sizeof(T));
    }
    if (buffer->length == buffer->capacity) {
        buffer->capacity *= 2;
        buffer->data = (T *) realloc(buffer->data, buffer->capacity * sizeof(T));
    }
    buffer->data[buffer->length++] = element;
}

template <typename T> void sbuff_free(StretchyBuffer<T> *buffer) {
    buffer->length = 0;
    buffer->capacity = 0;
    free(buffer->data);
    buffer->data = NULL;
}

void mesh_transform_vertex_positions_to_double_unit_box(int num_vertices, vec3 *vertex_positions) {
    vec3 L = V3(HUGE, HUGE, HUGE);
    vec3 R = V3(-HUGE, -HUGE, -HUGE);
    for (int i = 0; i < num_vertices; ++i) {
        L = cwiseMin(L, vertex_positions[i]);
        R = cwiseMax(R, vertex_positions[i]);
    }
    vec3 center = .5 * (L + R);
    vec3 size = R - L;
    double largest = MAX(MAX(size.x, size.y), size.z);
    for (int i = 0; i < num_vertices; ++i) {
        vertex_positions[i] -= center;
        vertex_positions[i] *= (2 / largest);
    }
}
BasicTriangleMesh3D load_basic_mesh(char *filename, bool transform_vertex_positions_to_double_unit_box) {
    BasicTriangleMesh3D basic_mesh = {};
    StretchyBuffer<vec3> vertex_positions = {};
    { // () load_basic_mesh
        FILE *fp = fopen(filename, "r");
        ASSERT(fp);
        char buffer[4096];
        while (fgets(buffer, NELEMS(buffer), fp) != NULL) {
            double x, y, z;
            ASSERT(sscanf(buffer, "%lf %lf %lf", &x, &y, &z) == 3);
            sbuff_push_back(&vertex_positions, { x, y, z });
        }
        fclose(fp);
    }
    basic_mesh.num_vertices = vertex_positions.length;
    basic_mesh.vertex_positions = vertex_positions.data; // NOTE stealing data pointer

    if (transform_vertex_positions_to_double_unit_box) {
        mesh_transform_vertex_positions_to_double_unit_box(basic_mesh.num_vertices, basic_mesh.vertex_positions);
    }
    return basic_mesh;
}

void fancy_mesh_alloc_compute_and_store_area_weighted_vertex_normals(FancyTriangleMesh3D *fancy_mesh) {
    ASSERT(fancy_mesh->vertex_normals == NULL);
    if (1) { // () fancy_mesh_alloc_compute_and_store_area_weighted_vertex_normals
        // TODO allocate fancy_mesh->vertex_normals        
        // TODO write entries of fancy_mesh->vertex_normals
        fancy_mesh->vertex_normals = (vec3 *) calloc(fancy_mesh->num_vertices, sizeof(vec3));
        for (int i_triangle = 0; i_triangle < fancy_mesh->num_triangles; ++i_triangle) {
            int3 ijk = fancy_mesh->triangle_indices[i_triangle];
            double A;
            vec3 n_hat;
            {
                vec3 abc[3];
                for (int d = 0; d < 3; ++d) {
                    abc[d] = fancy_mesh->vertex_positions[ijk[d]];
                }
                vec3 n = cross(abc[1] - abc[0], abc[2] - abc[0]);
                double mag_n = norm(n);
                A = norm(n) / 2;
                n_hat = n / mag_n;
            }
            for (int d = 0; d < 3; ++d) {
                fancy_mesh->vertex_normals[ijk[d]] += A * n_hat;
            }
        }
        for (int i_vertex = 0; i_vertex < fancy_mesh->num_vertices; ++i_vertex) {
            fancy_mesh->vertex_normals[i_vertex] = normalized(fancy_mesh->vertex_normals[i_vertex]);
        }
    }
}

FancyTriangleMesh3D load_fancy_mesh(char *filename, bool transform_vertex_positions_to_double_unit_box, bool compute_normals) {
    FancyTriangleMesh3D fancy_mesh = {};
    {
        StretchyBuffer<vec3> vertex_positions = {};
        StretchyBuffer<int3> triangle_indices = {};
        {
            FILE *fp = fopen(filename, "r");
            ASSERT(fp);
            char buffer[4096];
            while (fgets(buffer, NELEMS(buffer), fp) != NULL) {
                char prefix[16] = {};
                sscanf(buffer, "%s", prefix);
                if (strcmp(prefix, "f") == 0) {
                    int i, j, k;
                    ASSERT(sscanf(buffer, "%s %d %d %d", prefix, &i, &j, &k) == 4);
                    sbuff_push_back(&triangle_indices, { i - 1, j - 1, k - 1 });
                }
                if (strcmp(prefix, "v") == 0) {
                    double x, y, z;
                    ASSERT(sscanf(buffer, "%s %lf %lf %lf", prefix, &x, &y, &z) == 4);
                    sbuff_push_back(&vertex_positions, { x, y, z });
                }
            }
            fclose(fp);
        }
        // note: don't free the data pointers! (we're stealing them)
        fancy_mesh.num_triangles = triangle_indices.length;
        fancy_mesh.triangle_indices = triangle_indices.data;
        fancy_mesh.num_vertices = vertex_positions.length;
        fancy_mesh.vertex_positions = vertex_positions.data;
    }
    if (transform_vertex_positions_to_double_unit_box) {
        mesh_transform_vertex_positions_to_double_unit_box(fancy_mesh.num_vertices, fancy_mesh.vertex_positions);
    }
    if (compute_normals) {
        fancy_mesh_alloc_compute_and_store_area_weighted_vertex_normals(&fancy_mesh);
    }
    return fancy_mesh;
}

struct FPSCamera {
    vec3 origin;
    double angle_of_view;
    double theta;
    double phi;
};

mat4 fps_camera_get_C(FPSCamera *human) {
    return Translation(human->origin) * RotationY(human->theta) * RotationX(human->phi);
}

void fps_camera_move(FPSCamera *human) {
    vec3 ds = {}; {
        if (input.key_held['w']) { ds += transformVector(RotationY(human->theta), V3(0, 0, -1)); }
        if (input.key_held['s']) { ds += transformVector(RotationY(human->theta), V3( 0, 0, 1)); }
        if (input.key_held['a']) { ds += transformVector(RotationY(human->theta), V3(-1, 0, 0)); }
        if (input.key_held['d']) { ds += transformVector(RotationY(human->theta), V3( 1, 0, 0)); }
    }
    double norm_ds = norm(ds);
    if (!IS_ZERO(norm_ds)) {
        ds /= norm_ds;
        human->origin += ds;
    }
    if (glfwGetInputMode(window, GLFW_CURSOR) == GLFW_CURSOR_DISABLED) { // pointer lock
        human->theta -= input._mouse_dx_NDC;
        human->phi += input._mouse_dy_NDC;
        human->phi = CLAMP(human->phi, RAD(-80), RAD(80));
    }
}

void line_line_closest_points(vec3 a1, vec3 a2, vec3 b1, vec3 b2, vec3 *out_a_star, vec3 *out_b_star) {
    // http://www.geomalgorithms.com/algorithms.html#dist3D_Segment_to_Segment()
    // https://stackoverflow.com/questions/66979936/closest-two-3d-point-between-two-line-segment-of-varied-magnitude-in-different-p
    vec3 u = a2 - a1;
    vec3 v = b2 - b1;
    vec3 w = a1 - b1;
    double a = dot(u, u);         // always >= 0
    double b = dot(u, v);
    double c = dot(v, v);         // always >= 0
    double d = dot(u, w);
    double e = dot(v, w);
    double sc, sN, sD = a*c - b*b;  // sc = sN / sD, sD >= 0
    double tc, tN, tD = a*c - b*b;  // tc = tN / tD, tD >= 0
    double tol = 1e-15;
    // compute the line parameters of the two closest points
    if (sD < tol) {            // the lines are almost parallel
        sN = 0.0;              // force using point a1 on segment AB
        sD = 1.0;              // to prevent possible division by 0.0 later
        tN = e;
        tD = c;
    }
    else {                     // get the closest points on the infinite lines
        sN = (b*e - c*d);
        tN = (a*e - b*d);
    }
    // finally do the division to get sc and tc
    sc = (fabs(sN) < tol ? 0.0 : sN / sD);
    tc = (fabs(tN) < tol ? 0.0 : tN / tD);
    if (out_a_star) *out_a_star = a1 + (sc * u);
    if (out_b_star) *out_b_star = b1 + (tc * v);
}

void jank_widget_translate3D(mat4 PV, int num_points, vec3 *points) {
    if (widget_active_widget_ID != 0 && widget_active_widget_ID != WIDGET_ID_TRANSLATE) return;
    // please ignore; this function is jank
    static vec3 *selected_point;
    static vec3 *selected_handle;
    vec3 handle_colors[] = { monokai.white, monokai.white, monokai.white };
    double _L_handle_NDC = .1;
    double tol = .05;

    vec3 *hot = 0; {
        for (int i = 0; i < num_points; ++i) {
            vec3 *point = points + i;
            if (norm(transformPoint(PV, *point).xy - V2(input._mouse_x_NDC, input._mouse_y_NDC)) < tol) {
                hot = point;
            }
        }
    }

    static bool STILL_HOLDING_MOUSE_AFTER_SELECTING_POINT;
    if (!input.mouse_left_held) STILL_HOLDING_MOUSE_AFTER_SELECTING_POINT = false;
    if (hot && input.mouse_left_pressed) {
        STILL_HOLDING_MOUSE_AFTER_SELECTING_POINT = true;
        selected_point = (selected_point != hot) ? hot : 0;
    }

    if (hot) {
        basic_draw(POINTS, PV, 1, hot, !selected_point ? monokai.white : monokai.white, 10, true);
    }

    if (selected_point) {
        double L_handle = norm(*selected_point - transformPoint(inverse(PV), transformPoint(PV, *selected_point) + V3(_L_handle_NDC, 0, 0)));

        vec3 *hot_handle = 0;
        vec3 handles[3] = { *selected_point + V3(L_handle, 0, 0), *selected_point + V3(0, L_handle, 0), *selected_point + V3(0, 0, L_handle) };
        vec3 vertex_positions[] = { *selected_point, handles[0], *selected_point, handles[1], *selected_point, handles[2] };
        vec3 vertex_colors[] = { handle_colors[0], handle_colors[0], handle_colors[1], handle_colors[1], handle_colors[2], handle_colors[2] };
        basic_draw(LINES, PV, 6, vertex_positions, vertex_colors, 0);
        if (!STILL_HOLDING_MOUSE_AFTER_SELECTING_POINT) {
            for (int d = 0; d < 3; ++d) {
                if (norm(transformPoint(PV, handles[d]).xy - V2(input._mouse_x_NDC, input._mouse_y_NDC)) < tol) {
                    hot_handle = handles + d;
                }
            }
        }
        if (!selected_handle) {
            if (hot_handle) {
                basic_draw(POINTS, PV, 1, hot_handle, monokai.white, 6, true);
                if (input.mouse_left_pressed) {
                    selected_handle = hot_handle;
                }
            }
        } else {
            basic_draw(POINTS, PV, 1, selected_handle, monokai.white, 10, true);
        }
        if (input.mouse_left_held && selected_handle) {
            mat4 World_from_NDC = inverse(PV);
            vec3 s = transformPoint(World_from_NDC, V3(V2(input._mouse_x_NDC, input._mouse_y_NDC).x, V2(input._mouse_x_NDC, input._mouse_y_NDC).y, -1));
            vec3 t = transformPoint(World_from_NDC, V3(V2(input._mouse_x_NDC, input._mouse_y_NDC).x, V2(input._mouse_x_NDC, input._mouse_y_NDC).y, 1));
            vec3 new_handle_position;
            line_line_closest_points(*selected_point, *selected_handle, s, t, &new_handle_position, 0);
            *selected_point += new_handle_position - *selected_handle;

            {
                vec3 tmp = 16 * normalized(new_handle_position - *selected_point);
                vec3 tmp2[2] = { *selected_point - tmp, *selected_point + tmp };
                basic_draw(LINE_STRIP, PV, 2, tmp2, handle_colors[int(selected_handle - handles)], 2);
            }
            widget_active_widget_ID = WIDGET_ID_TRANSLATE;
        } else if (input.mouse_left_released) {
            selected_handle = 0;
            widget_active_widget_ID = 0;
        }
    }
}
