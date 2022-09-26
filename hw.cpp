#define COW_CRASH_ON_FLOATING_POINT_EXCEPTIONS
#define _CRT_SECURE_NO_WARNINGS
#include "snail.cpp"
#include "cow.cpp"
#include "_stretchy_buffer.cpp" // bare-bones templatized stretchy buffer

////////////////////////////////////////////////////////////////////////////////
// documentation                                                              //
////////////////////////////////////////////////////////////////////////////////

// mouse wheel to zoom

// HUGE is a real big number

// vec3 cross(vec3 a, vec3 b); // a x b  
// double squaredNorm(vec3 v); // |v|^2  
// double norm(vec3 v);        // |v|    
// vec3 normalized(vec3 v);    // v / |v|

// int3 is three contiguous ints  
//      and works a lot like vec3 
#if 0
int3 triangle = {};   // (0, 0, 0)
triangle.i = 4;       // (4, 0, 0)
triangle[1] = 5;      // (4, 5, 0)
triangle.data[2] = 6; // (4, 5, 6)
#endif                            

// // soup mesh                
//                             
// struct BasicTriangleMesh3D {
//     int num_vertices;       
//     vec3 *vertex_positions; 
// };                          

// // indexed mesh                      
//                                      
// struct FancyTriangleMesh3D {         
//     int num_vertices;                
//     int num_triangles;               
//     vec3 *vertex_positions;          
//     int3 *triangle_indices;          
//     vec3 *vertex_normals;            
// };                                   

// https://cplusplus.com/reference/cstdio/sscanf/ 
// https://cplusplus.com/reference/cstring/memcpy/

////////////////////////////////////////////////////////////////////////////////
// hw                                                                         //
////////////////////////////////////////////////////////////////////////////////

// begin please ignore these lines
void mesh_transform_vertex_positions_to_double_unit_box(int num_vertices, vec3 *vertex_positions);
void fancy_mesh_alloc_compute_and_store_area_weighted_vertex_normals(FancyTriangleMesh3D *fancy_mesh);
void fancy_mesh_merge_duplicated_vertices(FancyTriangleMesh3D *fancy_mesh);
// end please ignore these lines


// this function is already complete; you are free to use it as reference for load_basic_mesh(...)
FancyTriangleMesh3D load_fancy_mesh(char *filename, bool transform_vertex_positions_to_double_unit_box, bool compute_normals, bool unify_duplicated_vertices) {
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
                if (strcmp(prefix, "v") == 0) {
                    double x, y, z;
                    ASSERT(sscanf(buffer, "%s %lf %lf %lf", prefix, &x, &y, &z) == 4);
                    sbuff_push_back(&vertex_positions, { x, y, z });
                } else if (strcmp(prefix, "f") == 0) {
                    int i, j, k;
                    ASSERT(sscanf(buffer, "%s %d %d %d", prefix, &i, &j, &k) == 4);
                    sbuff_push_back(&triangle_indices, { i - 1, j - 1, k - 1 });
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
    if (unify_duplicated_vertices) {
        fancy_mesh_merge_duplicated_vertices(&fancy_mesh);
    }
    if (compute_normals) {
        fancy_mesh_alloc_compute_and_store_area_weighted_vertex_normals(&fancy_mesh);
    }
    return fancy_mesh;
}

// begin submission                                                             

BasicTriangleMesh3D load_basic_mesh(char *filename, bool transform_vertex_positions_to_double_unit_box) {
    BasicTriangleMesh3D basic_mesh = {};
    StretchyBuffer<vec3> vertex_positions = {};
    { // () load_basic_mesh
        // TODO populate vertex_positions using sbuff_push_back      
        // HINT see load_fancy_mesh(...) (or follow your heart)      
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

BasicTriangleMesh3D fancy2basic(FancyTriangleMesh3D fancy_mesh) {
    BasicTriangleMesh3D basic_mesh = {};
    { // () fancy2basic
        // TODO set basic_mesh.num_vertices         
        // TODO allocate basic_mesh.vertex_positions
        // TODO write basic_mesh.vertex_positions   
        fancy_mesh = fancy_mesh; // this is just here to temporarily remove unreferenced formal parameter compiler warning; please delete

    }
    return basic_mesh;
}

void mesh_transform_vertex_positions_to_double_unit_box(int num_vertices, vec3 *vertex_positions) {
    { // () mesh_transform_vertex_positions_to_double_unit_box
        // TODO overwrite entries of vetex_positions
        num_vertices = num_vertices; // this is just here to temporarily remove unreferenced formal parameter compiler warning; please delete
        vertex_positions = vertex_positions; // this is just here to temporarily remove unreferenced formal parameter compiler warning; please delete

    }
}

void fancy_mesh_alloc_compute_and_store_area_weighted_vertex_normals(FancyTriangleMesh3D *fancy_mesh) {
    ASSERT(fancy_mesh->vertex_normals == NULL);
    { // () fancy_mesh_alloc_compute_and_store_area_weighted_vertex_normals
        // TODO allocate fancy_mesh->vertex_normals        
        // TODO write entries of fancy_mesh->vertex_normals

    }
}

void fancy_mesh_merge_duplicated_vertices(FancyTriangleMesh3D *fancy_mesh) {
    int new_num_vertices = 0;
    vec3 *new_vertex_positions = (vec3 *) calloc(fancy_mesh->num_vertices, sizeof(vec3)); // (more space than we'll need)
    { // [] fancy_mesh_merge_duplicated_vertices
        // TODO set new_num_vertices
        // TODO wrie entries of new_vertex_positions                   
        // TODO overwrite entries of fancy_mesh->triangle_indices with new triangle indices                       
        // NOTE it is OK if your implementation is slow (mine takes ~5 seconds to fix up the teapot in debug mode)
        // NOTE please don't worry about space efficiency at all                                                  

    }
    if (new_num_vertices) {
        fancy_mesh->num_vertices = new_num_vertices;
        free(fancy_mesh->vertex_positions);
        fancy_mesh->vertex_positions = new_vertex_positions;
    }
}

void hw3a() {
    init();

    // no workarounds allowed :)

    BasicTriangleMesh3D basic_box = load_basic_mesh("data_basic_box", true);
    FancyTriangleMesh3D fancy_bunny = load_fancy_mesh("data_fancy_bunny", true, true, false);
    BasicTriangleMesh3D basic_bunny = fancy2basic(fancy_bunny);
    FancyTriangleMesh3D fancy_teapot_with_seams = load_fancy_mesh("data_fancy_teapot_with_seams", true, true, false);
    FancyTriangleMesh3D fancy_teapot_no_seams = load_fancy_mesh("data_fancy_teapot_with_seams", true, true, true);

    Camera3D camera = { 5, RAD(45), RAD(0), RAD(-10), 0, .1 };
    double t = 0;
    bool paused = false;
    int part = 0;
    while (begin_frame()) {
        camera_move(&camera);
        mat4 P = camera_get_P(&camera);
        mat4 V = camera_get_V(&camera);
        mat4 M = RotationY(.7*t);
        mat4 PVM = P * V * M;

        int num_parts = 0;
        if (part == num_parts++) { basic_draw(TRIANGLE_MESH, PVM, basic_box, monokai.green, AVG(monokai.green, monokai.white), 3); }
        if (part == num_parts++) { fancy_draw(P, V, M, fancy_bunny, monokai.blue); }
        if (part == num_parts++) { basic_draw(TRIANGLE_MESH, PVM, basic_bunny, monokai.blue, AVG(monokai.blue, monokai.white), 3); }
        if (part == num_parts++) { fancy_draw(P, V, M, fancy_teapot_with_seams, monokai.red); }
        if (part == num_parts++) { fancy_draw(P, V, M, fancy_teapot_no_seams, monokai.red); }
        imgui_slider("part", &part, 0, num_parts - 1, 'j', 'k', true);

        { // bounding [-1, 1]^3 box
            double tmp[] = {-1,-1,-1,-1,1,-1,-1,1,1,-1,-1,1,-1,-1,-1,1,-1,-1,1,1,-1,1,1,1,1,-1,1,1,-1,-1,1,1,-1,-1,1,-1,1,1,-1,1,1,1,-1,1,1,1,1,1,1,-1,1,-1,-1,1,1,-1,1,};
            basic_draw(LINE_STRIP, PVM.data, XYZ, RGB, NELEMS(tmp) / 3, tmp, NULL, 1, 1, 1, .1, 2);
        }

        imgui_checkbox("paused", &paused, 'p');
        if (!paused) { t += .0167; }
        if (imgui_button("reset", 'r')) { t = 0; }
    }
}

void hw3b() {
    init();

    mat4 M[10] = {}; {
        for (int i = 0; i < NELEMS(M); ++i) {
            M[i] = Identity4x4;
        }
    }

    vec2 L_quads[] = {{0,0},{2,0},{2,1},{0,1},{0,0},{1,0},{1,3},{0,3}};

    Camera2D camera = { 20 };
    while (begin_frame()) {
        camera_move(&camera);
        mat4 PV = camera_get_PV(&camera);

        // // 2D transform API                                                                     
        // mat4 Translation(double t_x, double t_y); // translation by (t_x, t_y)^T                
        // mat4 Scaling(double s_x, double s_y) // scaling by s_x in x and x_y in y                
        // mat4 Rotation(double theta); // counter-clockwise rotation about origin by theta radians
        // - NOTE RAD(deg) converts from degrees to radians                                        

        #if 0
        M[0] = Translation(5, 5);
        #endif

        // () yellow
        M[0] = {};

        // () purple
        M[1] = {};

        // () orange
        M[2] = {};

        // () lightblue
        M[3] = {};

        // () red
        M[4] = {};

        // () buff (tan)
        M[5] = {};

        // () gray
        M[6] = {};

        // [] green
        M[7] = {};

        // [] purplishpink
        // NOTE only a perfect solution will score full credit                         
        //      i.e., no hard-coded constants with 10 digits after the decimal place :)
        M[8] = {};

        // <> blue
        // NOTE only a perfect solution will score full credit                         
        //      i.e., no hard-coded constants with 10 digits after the decimal place :)
        M[9] = {};

        { // draw L-block's
            basic_draw(QUADS, PV, NELEMS(L_quads), L_quads, monokai.white);
            for (int i = 0; i < NELEMS(M); ++i) {
                basic_draw(QUADS, PV * M[i], NELEMS(L_quads), L_quads, color_get_kelly(i));
            }
        }
        { // bespoke widget
            vec2 mouse_position = input_get_mouse_position_in_world_coordinates(PV);
            int x = (int) roundf((float) mouse_position.x);
            int y = (int) roundf((float) mouse_position.y);
            imgui_readout("x", &x);
            imgui_readout("y", &y);
            gl_PV(PV);
            gl_color(0, 1, 0, 1);
            gl_begin(GL_POINTS);
            gl_vertex(x, y);
            gl_end();
            gl_begin(GL_LINE_LOOP);
            gl_vertex(x, y);
            gl_color(0, 1, 0, .5);
            gl_vertex(0, y);
            gl_vertex(0, 0);
            gl_vertex(x, 0);
            gl_end();
        }

        // NOTE if you want to draw other stuff to help you debug, do it down here

    }
}

// // 3D transform API                                               
// mat4 Translation(double t_x, double t_y, double t_z);             
// mat4 Scaling(double s_x, double s_y, double s_z);                 
// mat4 RotationX(double theta);                                     
// mat4 RotationY(double theta);                                     
// mat4 RotationZ(double theta);                                     
// mat4 Rotation(vec3 axis, dbouel theta); // probably not so useful?

void hw3c() {
    init();
    Camera3D camera = { 10, RAD(45) };
    double t = 0;
    bool playing = false;
    while (begin_frame()) {                                    
        camera_move(&camera);
        mat4 P = camera_get_P(&camera);
        mat4 V = camera_get_V(&camera);

        mat4 S = Scaling(1.25 - .25 * cos(5 * t), .722 + .278 * cos(5 * t), 1.25 - .25 * cos(5 * t));
        fancy_draw(P, V, Translation(-4.5, 0, 0) * S, meshlib.fancy_box, monokai.red);
        fancy_draw(P, V, Translation(-1.5, 0, 0) * S, meshlib.fancy_cone, monokai.yellow);
        fancy_draw(P, V, Translation( 1.5, 0, 0) * S, meshlib.fancy_cylinder, monokai.blue);
        fancy_draw(P, V, Translation( 4.5, 0, 0) * S, meshlib.fancy_sphere, monokai.purple);

        { // floor
            gl_PV(P * V);
            gl_begin(QUADS);
            gl_color(1, 1, 1, .5);
            gl_vertex(10, 0, 10);
            gl_vertex(-10, 0, 10);
            gl_vertex(-10, 0, -10);
            gl_vertex(10, 0, -10);
            gl_end();
        }

        imgui_checkbox("playing", &playing, 'p');
        if (playing) { t += .0167; }
    }
}

void hw() {
    init(false, "", 540, 0, 100);
    hw3a();
    hw3b();
    hw3c();
}

// end submission                                                               


int main() {
    hw();
    return 0;
}

