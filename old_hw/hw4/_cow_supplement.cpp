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
