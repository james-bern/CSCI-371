#version 330 core

uniform float iTime;
uniform vec2 iResolution;
uniform mat4 C;

out vec4 fragColor;

float TAU = 6.28318530718;
float RAD(float deg) {
    return deg / 360.0 * TAU;
}

float sdTorus(vec3 p, vec2 t) {
    vec2 q = vec2(length(p.xz)-t.x,p.y);
    return length(q)-t.y;
}
float sdBox(vec3 p, vec3 b) {
    vec3 q = abs(p) - b;
    return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

vec4 march(vec3 o, vec3 d) {
    // params
    const int   march_max_steps     = 80;
    const float march_hit_tolerance = 0.001;
    const float march_max_distance  = 100.0;
    // o -- camera origin               
    // t -- distance marched along ray  
    // d -- camera direction            
    // p -- current position along ray  
    // f -- distance to implicit surface
    float t = 0.0;
    int step = 0;
    while ((step++ < march_max_steps) && (t < march_max_distance)) {
        vec3 p = o + t * d;
        float f = march_max_distance; {
            // TODO: make the scene of your dreams :)
            { // box
                vec3 box_position = vec3(2.0 * sin(iTime), 0.0, 0.0);
                vec3 box_side_lengths = vec3(1.0);
                float distance_to_box = sdBox(p - box_position, box_side_lengths);
                f = min(f, distance_to_box);
            }
            { // torus
                vec3 torus_position = vec3(0.0, 0.0, 0.0);
                float torus_major_radius = 2.0;
                float torus_minor_radius = 0.25;
                vec2 torus_radii = vec2(torus_major_radius, torus_minor_radius);
                float distance_to_torus = sdTorus(p - torus_position, torus_radii);
                f = min(f, distance_to_torus);
            }
        }
        if (f < march_hit_tolerance) { // hit!
            return vec4(0.5 + 0.5 * cos(TAU * (vec3(0.0, 0.33, -0.33) - vec3(0.3 * p.z))), 1.0);
        }
        t += min(f, .5); // make the number smaller if you're getting weird artifacts
    }
    return vec4(0.0);
}

void main() {
    vec3 o = C[3].xyz; // glsl is pretty neato :)
    vec3 d; {
        float theta_over_two = RAD(30.0);
        vec2 d_xy_camera = (gl_FragCoord.xy - (0.5 * iResolution.xy)) * (tan(theta_over_two) / (0.5 * iResolution.y));
        vec3 d_camera = normalize(vec3(d_xy_camera, -1.0));
        d = mat3(C) * d_camera;
    }
    fragColor = march(o, d);
}
