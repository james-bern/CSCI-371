// o refers to the link'o center of mass

//  0     1          3 --> x


//  ---o--J----o-----       
//   ^link0    ^link1       

struct Link {
    vec2 o;
    real theta;
    real L; // for vis only
};


// optimization vector
/*
   o0
   o1
   o2
   ...
   theta0
   theta1
   theta2
   ...
 */


const int num_links = 2;
double optimization_vector[] = {
    0.5,
    0.0,
    2.0,
    0.0,

    0.0,
    RAD(45),
};

vec2 *link_o;
real *link_theta;
double link_L[] = { 1.0, 2.0 };




struct Joint {
    int i[2];  // joint is connected to links[this->io]
    vec2 S[2]; // position in coordinate system of links[i0]
};
const int num_joints = 1;
Joint joints[num_joints] = {
    {
        { 0, 1 },
        { {  0.5, 0.0 }, { -1.0, 0.0 } },
    },
};

void app_linkage_simulation() {
    Camera2D camera = { 10.0 };

    link_o = (vec2 *) optimization_vector;
    link_theta = (real *) (optimization_vector + 2 * num_links);

    while (cow_begin_frame()) {

        // todo solve physics

        { // gui and draw
            mat4 PV = camera_get_PV(&camera);
            widget_drag(PV, num_links, (vec2 *) optimization_vector);
            for_(link_i, num_links) {
                static char buffer[16];
                sprintf(buffer, "theta_%d", link_i);
                gui_slider(buffer, &link_theta[link_i], -PI, PI, true);
                sprintf(buffer, "x_%d", link_i);
                gui_readout(buffer, &link_o[link_i].x);
                sprintf(buffer, "y_%d", link_i);
                gui_readout(buffer, &link_o[link_i].y);
            }

            eso_begin(PV, SOUP_LINES); {
                for_(link_i, num_links) {
                    eso_color(color_kelly(link_i));
                    for_sign(sign) {
                        eso_vertex(link_o[link_i] + sign * .5 * link_L[link_i] * e_theta(link_theta[link_i]));
                    }
                }
            } eso_end();

            eso_begin(PV, SOUP_POINTS); {

                eso_color(monokai.blue);

                for_(joint_i, num_joints) {
                    Joint *joint = joints + joint_i;
                    for_(d, 2) {

                        vec2 S = joint->S[d];

                        int link_i = joint->i[d];
                        vec2 o = link_o[link_i];
                        vec2 x = e_theta(link_theta[link_i]);
                        vec2 y = e_theta(link_theta[link_i] + RAD(90));


                        eso_vertex(o + S.x * x + S.y * y);

                    } 
                }
            } eso_end();
        }
    }
}








