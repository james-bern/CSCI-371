void app_space_fish_arduino() {
    // turn on/off led from cow
    // turn on/off virtual led Arduino
    // closed loop mass-spring control to hit target axle_theta
    while (cow_begin_frame()) {

    }
}

void app_space_fish_2D() {
    vec2 s_motor = {};
    vec2 s_magnet = s_motor + V2(100.0, 0.0);
    real theta_motor = RAD(90);
    real theta_magnet = RAD(90);

    Camera2D camera = { 250.0 };
    double timestep = .0167;
    double time = -timestep;
    bool playing = false;
    while (cow_begin_frame()) {
        gui_checkbox("playing", &playing, 'p');
        time += timestep;
        camera_move(&camera);
        mat4 PV = camera_get_PV(&camera);

        if (playing) {
            theta_motor += RAD(sin(3.0 * time));
            theta_magnet += RAD(sin(3.0 * time));
        }
        gui_slider("theta_motor", &theta_motor, RAD(0), RAD(270));
        gui_slider("theta_magnet", &theta_magnet, RAD(0), RAD(270));

        #define FISH_STATE_NORMAL_TRACKING 0
        #define FISH_STATE_PUSHED 1
        #define FISH_STATE_ERROR 2
        int fish_state = FISH_STATE_NORMAL_TRACKING;
        double Delta = theta_magnet - theta_motor;
        double r_tol = RAD(5);
        if (Delta > r_tol) {
            fish_state = FISH_STATE_PUSHED;
        } else if (Delta < -r_tol) {
            fish_state = FISH_STATE_ERROR;
        }

        vec2 s_third = s_motor + 50.0 * e_theta(theta_motor);
        vec2 s_fourth = s_magnet + 50.0 * e_theta(theta_magnet);
        vec2 s_fish = s_magnet + 100.0 * e_theta(theta_magnet);
        eso_begin(PV, SOUP_LINES); {
            eso_color(monokai.blue);
            eso_vertex(s_motor);
            eso_vertex(s_magnet);
            eso_color(monokai.orange);
            eso_vertex(s_motor);
            eso_vertex(s_third);
            eso_vertex(s_magnet);
            eso_vertex(s_fish);
            eso_color(monokai.purple);
            eso_vertex(s_third);
            eso_vertex(s_fourth);
        } eso_end();
        eso_begin(PV, SOUP_LINES); {
            eso_color(
                    (fish_state == FISH_STATE_NORMAL_TRACKING) ? monokai.orange :
                    (fish_state == FISH_STATE_PUSHED) ? monokai.blue :
                    monokai.red);
            int N = 64;
            for_(i, N) {
                eso_vertex(s_fish + 25 * e_theta(NUM_DENm1(i, N) * 2.0 * PI));
            }
        } eso_end();
    }
}

void app_sapce_fish_exploration() {
    real dynamixel_Z_theta = 0.0;
    real dynamixel_Y_theta = 0.0;
    real asm5601_theta = 0.0;

    real axle_theta = 0.0;

    // todo save camera to camera.txt super easily
    Camera3D camera = { 10.0 };
    while (cow_begin_frame()) {
        camera_move(&camera);
        mat4 P, V, PV; {
            P = camera_get_P(&camera);
            V = camera_get_V(&camera);
            PV = P * V;
        }

        gui_slider("dynamixel_Z_theta", &dynamixel_Z_theta, RAD(-180), RAD(180), true);
        gui_slider("dynamixel_Y_theta", &dynamixel_Y_theta, RAD(-180), RAD(180), true);
        gui_slider("asm5601_theta", &asm5601_theta, RAD(-180), RAD(180), true);

        int kelly_i = 0;
        auto _Q = [&](vec3 size, vec3 origin_datum, mat4 M2) {
            // // origin_datum
            // (0, 0, 0) center
            // (0, -1, 0) center of bottom face
            // (1, 1, 0) center of upper right edge
            // (1, 1, 1) ... corner
            mat4 M = M4_Translation(-cwiseProduct(origin_datum, .5 * size)) * M4_Scaling(.5 * size);
            library.meshes.box.draw(P, V, M2 * M, AVG(monokai.gray, color_kelly(kelly_i++)));
            // library.soups.box.draw(P * V * M2 * M, color_kelly(kelly_i++));
        };
        auto Q = [&](real size_x, real size_y, real size_z, real datum_x, real datum_y, real datum_z, mat4 M2 = globals.Identity) {
            _Q({ size_x, size_y, size_z }, { datum_x, datum_y, datum_z }, M2);
        };

        axle_theta += .0167;
        vec3 s_virtual = { 4.0 * sin(axle_theta), 4.0 };
        dynamixel_Y_theta = LERP(CLAMP(INVERSE_LERP(s_virtual.x, -1.0, 1.0), 0.0, 1.0), -PI, 0);

        vec3 S_axle = { 3.0, 1.5, 0.0 };
        vec3 s_axle = transformPoint(M4_RotationAboutYAxis(dynamixel_Y_theta), S_axle);


        // todo transform heirarchy

        double target_asm5601_theta = atan2((s_virtual - s_axle).xy);
        asm5601_theta = CLAMP(target_asm5601_theta, RAD(60), RAD(120));

        real L_bar = 2.5;
        vec3 s_real = s_axle + L_bar * V3(e_theta(asm5601_theta), 0.0);


        {
            eso_begin(PV, SOUP_LINES, 1.0, true);
            eso_color(monokai.blue);
            eso_vertex(s_axle);
            eso_vertex(s_real);
            eso_end();
        }
        library.meshes.sphere.draw(P, V, M4_Translation(s_axle) * M4_Scaling(0.1), monokai.blue);
        library.meshes.sphere.draw(P, V, M4_Translation(s_real), monokai.blue);
        library.meshes.sphere.draw(P, V, M4_Translation(s_virtual), monokai.orange);


        Q(10.0, 0.1, 10.0, 0.0, -1.0, 0.0, M4_RotationAboutYAxis(dynamixel_Y_theta));
        // Q(1.0, 3.0, 0.1, 0.0, -1.0, 0.0, M4_Translation(1.0, 1.0) * M4_RotationAboutZAxis(-dynamixel_Z_theta));
        // Q(1.0, 3.0, 0.1, 0.0, -1.0, 0.0, M4_Translation(3.0, 1.0) * M4_RotationAboutZAxis(-asm5601_theta));


    }
}
