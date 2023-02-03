// the artist prescribes a target motion
// the system displays the 

void app_space_fish_2D() {
    #define ACTUAL_FISH_STATE_NORMAL 0
    #define ACTUAL_FISH_STATE_PUSHED 1
    #define ACTUAL_FISH_STATE_ERROR 2
    #define ACTUAL_FISH_STATE_HIDDEN 3
    vec3 actual_fish_state_color[] = { monokai.orange, monokai.blue, monokai.red, monokai.green };

    vec2 s_motor = {};
    vec2 s_magnet = s_motor + V2(100.0, 0.0);
    real theta_motor = RAD(90);
    real theta_magnet = RAD(90);

    real r_fish = 25.0;

    real x_cutoff = 80.0;

    vec2 s_virtual = {};

    vec2 s_dragger = { 50.0, 100.0 };
    vec2 s_target;

    // ik is always fun
    // maybe the animator is specifying target behaviors
    // - a chance to work with someone in dance could be fun
    // - lion tamer?
    // - bullfight

    // i think we actually do want a simulation
    // funny

    // well let's just model it as a mass spring system
    // do we end up back at impedance control?

    // what can we do without a simulation?

    // what assumptions can we make?
    // - human is infinitely stiff and strong

    // taking a step back questions
    // - do we want/need the silicone plate?
    // - if not, then the human can touch the fish any time any where
    // - - this seems like a very different (albeit interesting) problem
    // - - also not a fish


    Camera2D camera = { 250.0 };
    double timestep = .0167;
    double time = -timestep;
    bool playing = false;
    while (cow_begin_frame()) {
        camera_move(&camera);
        mat4 PV = camera_get_PV(&camera);

        gui_slider("theta_motor", &theta_motor, RAD(0), RAD(270));
        gui_slider("theta_magnet", &theta_magnet, RAD(0), RAD(270));
        gui_checkbox("playing", &playing, 'p');
        if (playing) {
            time += timestep;
            theta_magnet = RAD(90.0 + 10.0 * sin(3.0 * time));
        }

        double Delta = theta_magnet - theta_motor;
        double r_tol = RAD(5);

        vec2 s_actual = s_magnet + 100.0 * e_theta(theta_magnet);
        int fish_state; {
            if (s_actual.x + r_fish < x_cutoff) {
                fish_state = ACTUAL_FISH_STATE_HIDDEN;
            } else {
                fish_state = ACTUAL_FISH_STATE_NORMAL;
                if (Delta > r_tol) {
                    fish_state = ACTUAL_FISH_STATE_PUSHED;
                } else if (Delta < -r_tol) {
                    fish_state = ACTUAL_FISH_STATE_ERROR;
                }
            }
        }

        if (fish_state != ACTUAL_FISH_STATE_HIDDEN) {
            s_target = s_virtual = s_actual; // *
        } else {
            s_target = s_dragger;
            s_virtual += magClamped(s_target - s_virtual, 3.0);
        }

        { // draw
            eso_begin(PV, SOUP_LINES); {
                vec2 _s_third = s_motor + 50.0 * e_theta(theta_motor);
                vec2 _s_fourth = s_magnet + 50.0 * e_theta(theta_magnet);
                eso_color(monokai.gray);
                eso_vertex(s_motor);
                eso_vertex(s_magnet);
                eso_vertex(s_motor);
                eso_vertex(_s_third);
                eso_vertex(s_magnet);
                eso_vertex(s_actual);
                eso_color(monokai.purple);
                eso_vertex(_s_third);
                eso_vertex(_s_fourth);
            } eso_end();
            int N = 64;
            eso_begin(PV, SOUP_TRIANGLE_FAN); {
                eso_color(AVG(monokai.white, actual_fish_state_color[fish_state]));
                for_(i, N) {
                    eso_vertex(s_virtual + r_fish * e_theta(NUM_DENm1(i, N) * 2.0 * PI));
                }
            } eso_end();
            eso_begin(PV, SOUP_LINE_LOOP); {
                eso_color(actual_fish_state_color[fish_state]);
                for_(i, N) {
                    eso_vertex(s_actual + r_fish * e_theta(NUM_DENm1(i, N) * 2.0 * PI));
                }
            } eso_end();
            eso_begin(PV, SOUP_POINTS, 24.0); {
                eso_color(actual_fish_state_color[fish_state]);
                eso_vertex(s_actual + V2(r_fish, 0));
            } eso_end();
            eso_begin(PV, SOUP_LINES); {
                eso_color(monokai.white, 0.5);
                double r = 15.0;
                eso_vertex(x_cutoff, 100.0 + r);
                eso_vertex(x_cutoff, 100.0 - r);
            } eso_end();
        }

        widget_drag(PV, 1, &s_dragger);
        soup_draw(PV, SOUP_POINTS, 1, &s_dragger, NULL, V3(1.0, 0.0, 1.0));
        soup_draw(PV, SOUP_POINTS, 1, &s_target, NULL, V3(1.0, 0.0, 1.0), 24.0);
    }
}

void app_space_fish_electronics() {
    // turn on/off led from cow
    // turn on/off virtual led from Arduino
    // bluetooth
    // dynamixel
    // magnet
    // closed loop model-based control to hit target axle_theta
    while (cow_begin_frame()) {

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
