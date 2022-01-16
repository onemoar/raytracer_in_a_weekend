#define SPHERES_ONLY 0
static int samples_per_pixel = 2500;

struct ray {
    vec3 origin;
    vec3 direction;

    vec3 at(float t) {
        return origin + t * direction;
    }
};


struct camera {
    camera() {}
    camera(vec3 pos, vec3 direction, vec3 world_up, float fov, float aspect_ratio, float aperture, float focus_dist) {
        float width = tan(fov * DEG_TO_RAD / 2);
        float viewport_width = 2 * width;
        float viewport_height = viewport_width / aspect_ratio;

        forward = normalize(-direction);
        right = normalize(cross(world_up, forward));
        up = cross(forward, right);

        position = pos;

        horizontal = focus_dist * viewport_width * right;
        vertical   = focus_dist * viewport_height * up;
        lower_left_corner = position - horizontal/2 - vertical/2 - forward;

        lens_radius = aperture / 2;
    }

    ray get_ray(float u, float v) {
        vec3 rd = lens_radius * random_in_unit_disk();
        vec3 offset = right * rd.x + up * rd.y;

        return ray {position + offset, normalize(lower_left_corner + u * horizontal + v * vertical - position - offset)};
    }

    vec3 position;
    vec3 horizontal;
    vec3 vertical;
    vec3 lower_left_corner;
    vec3 forward;
    vec3 right;
    vec3 up;

    float lens_radius;
};

struct material;

struct ray_hit {
    vec3 position;
    vec3 normal;
    double t;
    material* material;
    bool front_face;

    inline void set_face_normal(ray& r, vec3& outward_normal) {
        front_face = dot(r.direction, outward_normal) < 0;
        normal = front_face ? outward_normal :-outward_normal;
    }
};

struct sphere {
    vec3 position;
    float radius;

    material* material;
};

struct AABB {
    vec3 min;
    vec3 max;

    material* material;
};

struct cone {
    vec3 position;
    vec3 direction;
    float height;
    float cos_angle_sq;

    material* material;
};

struct disc {
    vec3 position;
    vec3 normal;
    float radius;

    material* material;
};

camera current_camera;
std::vector<sphere> spheres;
std::vector<AABB> boxes;
std::vector<cone> cones;
std::vector<disc> discs;

enum struct material_type {
    diffuse,
    metal,
    dielectric
};

struct material {
    vec3 color;
    float roughness;
    float refraction_index;
    material_type type;
    bool is_emissive = false;

    material() {}
    material(vec3 c, float r, float ri, material_type t, bool e) {
        color = c;
        roughness = r;
        refraction_index = ri;
        type = t;
        is_emissive = e;
    }
};

float reflectance(float cosine, float ri) {
    // Use Schlick's approximation for reflectance.
    auto r0 = (1-ri) / (1+ri);
    r0 = r0*r0;
    return r0 + (1-r0)*pow((1 - cosine),5);
}

bool scatter(ray& in, ray_hit& hit, vec3& c, ray& scattered) {
    material mat = *hit.material;

    switch(mat.type) {
        case material_type::diffuse: {
            c = mat.color;

            if (mat.is_emissive) return false;

            vec3 scatter_dir = hit.normal + random_on_unit_sphere();

            if (length(scatter_dir) < 1e-7)
                scatter_dir = hit.normal;

            scattered = ray {hit.position, normalize(scatter_dir)};

            return true;
        } break;

        case material_type::metal: {
            vec3 reflected = reflect(normalize(in.direction), hit.normal);

            scattered = ray {hit.position, normalize(reflected + mat.roughness * random_on_unit_sphere())};
            c = mat.color;

            return true;
        } break;

        case material_type::dielectric: {
            c = mat.color;
            double refraction_ratio = hit.front_face ? (1.0/mat.refraction_index) : mat.refraction_index;

            vec3 unit_direction = normalize(in.direction);
            double cos_theta = fmin(dot(-unit_direction, hit.normal), 1.0);
            double sin_theta = sqrt(1.0 - cos_theta*cos_theta);

            bool cannot_refract = refraction_ratio * sin_theta > 1.0;
            vec3 direction;
            if (cannot_refract || reflectance(cos_theta, refraction_ratio) > random_float_01())
                direction = reflect(unit_direction, hit.normal);
            else
                direction = refract(unit_direction, hit.normal, refraction_ratio);

            scattered = ray {hit.position, normalize(direction)};
            return true;
        } break;
    }

    return false;
}

double hit_sphere(sphere& s, ray& r, float& a) {
    vec3 oc = r.origin - s.position;

    auto half_b = dot(oc, r.direction);
    auto c = dot(oc, oc) - s.radius*s.radius;

    auto discriminant = half_b*half_b - c;
    if (discriminant < 0) {
        return -1.0;
    } else {
        return (-half_b - sqrt(discriminant));
    }
}

double hit_cone(cone& cone, ray& r) {
    auto oc = r.origin - cone.position;
    auto bd = dot(r.direction, cone.direction);
    auto ad = dot(oc, cone.direction);
    auto dd = dot(cone.direction, cone.direction) * cone.cos_angle_sq;
    auto ab = dot(oc, r.direction);
    auto aa = dot(oc, oc);

    auto a = bd * bd - dd;
    // auto a = bd * bd - dd;
    auto half_b = ad * bd - ab * dd;
    auto c = ad * ad - aa * dd;

    auto discriminant = half_b*half_b - a*c;
    if (discriminant < 0) {
        return -1.0;
    } else {
        return (-half_b + sqrt(discriminant)) / a;  // почему + а не - ???
    }
}

double hit_aabb(AABB& a, ray& r, vec3& normal) {
    double tmin = 0;
    double tmax = 9999999999;
    for (int i = 0; i < 3; i++) {
        if (fabs(r.direction[i]) < 1e-7) {
            if (r.origin[i] < a.min[i] || r.origin[i] > a.max[i]) return -1.0;
        } else {
            float ood = 1.0f / r.direction[i];
            float t1 = (a.min[i] - r.origin[i]) * ood;
            float t2 = (a.max[i] - r.origin[i]) * ood;

            bool didSwap = t1 > t2;
            if (didSwap) 
                swap(t1, t2);
            if (t1 > tmin) {
                tmin = t1;

                // TODO: better way to calcuclate normal?
                normal = v3_zero;
                normal[i] = (1 - 2 * (int)didSwap);
            }
            if (t2 < tmax) tmax = t2;

            if (tmin > tmax) return -1.0;
        }
    }

    return tmin;
}

double hit_disc(disc& d, ray& r) {
    auto nb = dot(d.normal, r.direction);
    if (fabs(nb) < 1e-3)
        return -1;
    auto np = dot(d.normal, d.position);
    auto na = dot(d.normal, r.origin);

    float t = (np - na) / nb;
    auto distance = length(r.at(t) - d.position);
    if (distance > d.radius) return -1;

    return t;
}

uint64 ray_count = 0;
bool raycast(ray& r, ray_hit& hit) {
    ray_count++;

    double closest_t = 9999999999;
    int sphere_id = -1;
    float epsilon = 1e-3;
    auto a = dot(r.direction, r.direction);

    for (int i = 0; i < spheres.size(); i++) {
        auto t = hit_sphere(spheres[i], r, a);
        if (t > epsilon && t < closest_t) {
            closest_t = t;
            sphere_id = i;
        }
    }

    #if !SPHERES_ONLY
    int box_id = -1;

    vec3 box_normal;
    for (int i = 0; i < boxes.size(); i++) {
        auto t = hit_aabb(boxes[i], r, box_normal);
        if (t > epsilon && t < closest_t) {
            closest_t = t;
            box_id = i;
        }
    }

    int cone_id = -1;
    for (int i = 0; i < cones.size(); i++) {
        auto t = hit_cone(cones[i], r);
        if (t > epsilon && t < closest_t) {
            auto cone = cones[i];
            auto proj = dot(r.at(t) - cone.position, cone.direction);
            if (proj > 0 && proj < cone.height) {
                closest_t = t;
                cone_id = i;
            }
        }
    }

    int disc_id = -1;
    for (int i = 0; i < discs.size(); i++) {
        auto t = hit_disc(discs[i], r);
        if (t > epsilon && t < closest_t) {
            closest_t = t;
            disc_id = i;
        }
    }

    if (disc_id != -1) {
        hit.t = closest_t;
        hit.position = r.at(closest_t);
        hit.set_face_normal(r, discs[disc_id].normal);
        hit.material = discs[disc_id].material;
        return true;
    }
    if (cone_id != -1) {
        hit.t = closest_t;
        hit.position = r.at(closest_t);

        //@TODO: this is wrong and non diffuse materials will not work properly
        hit.normal = random_on_unit_sphere();
        hit.material = cones[cone_id].material;        
        return true;
    }
    if (box_id != -1) {
        hit.t = closest_t;
        hit.position = r.at(closest_t);
        hit.set_face_normal(r, box_normal);
        hit.material = boxes[box_id].material;        
        return true;
    }

    #endif
    if (sphere_id != -1) {
        hit.t = closest_t;
        hit.position = r.at(closest_t);
        vec3 normal = (hit.position - spheres[sphere_id].position) / spheres[sphere_id].radius;
        hit.set_face_normal(r, normal);
        hit.material = spheres[sphere_id].material;
        return true;
    }

    return false;
}

vec3 get_pixel_color(camera& camera, int x, int y) {
    vec3 result = {};

    for (int i = 0; i < samples_per_pixel; i++) {
        float u = (x + random_float_01()) * 1.0 / (bitmap_width - 1);            
        float v = (bitmap_height - (y + random_float_01())) * 1.0 / (bitmap_height - 1);            

        ray r = camera.get_ray(u, v);

        vec3 color_acumulator = v3_one;
        for (int j = 0; j < 50; j++) {
            ray_hit hit;
            if (raycast(r, hit)) {
                vec3 hit_color;
                ray scattered;

                if (scatter(r, hit, hit_color, scattered)) {
                    color_acumulator = color_acumulator * hit_color;
                    r = scattered;
                } 
                else if (hit.material->is_emissive) {
                    color_acumulator = color_acumulator * hit_color;
                    break;
                }
            } else {
                // vec3 dir = normalize(r.direction);
                // float sky_t = 0.5 * (dir.y + 1);
                // color_acumulator = color_acumulator * ((1 - sky_t) * vec3(1,1,1) + sky_t * vec3(.5, .7, 1)); 
                color_acumulator = v3_zero;
                break;
            }
        }
        result += color_acumulator;
    }

    result /= samples_per_pixel;
    result.x = sqrt(result.x);
    result.y = sqrt(result.y);
    result.z = sqrt(result.z);
    if (result.x > 1) result.x = 1;
    if (result.y > 1) result.y = 1;
    if (result.z > 1) result.z = 1;

    return result;
}

static bool render_scene_row(camera& camera, int current_row) {
    int pitch = bitmap_width*bytes_per_pixel;
    uint8 *row_start = (uint8 *)bitmap_buffer + current_row * pitch;    

    for(int y = current_row; y < bitmap_height; y++) {
        uint32 *pixel = (uint32 *)row_start;
        for(int x = 0; x < bitmap_width; x++) {
            vec3 c = get_pixel_color(camera, x, y);

            uint8 red = (int)(c.r * 255);
            uint8 green = (int)(c.g * 255);
            uint8 blue = (int)(c.b * 255);

            *pixel++ = (red << 16) | (green << 8) | blue;
        }

        row_start += pitch;
        return false;
    }

    return true;
}

void add_spiraling_spheres(vec3 position, float angle, float height, int num_balls) {
    float current_height = 0;
    for (int i = 0; i < num_balls; i++) {
        float radius = current_height * angle;
        current_height = - i * height / num_balls;
        auto random = random_float(0, 1);
        material* mat;
        auto ball_radius = random_float(.05,.1);
        if (random < .5) {
            int channel = xor_shift_32() % 2;
            vec3 color = {0,0,0};
            color[channel] = 2;
            mat = new material(color, 0, 0, material_type::diffuse, true);
            ball_radius = .05;
        } else 
        if (random < .75) {
            mat = new material(vec3_random(0,1), random_float(0,0.1), 0, material_type::metal, false);
        } else {
            mat = new material(vec3_random(0,1), 0, 1 + random_float(0,1), material_type::dielectric, false);
        }
        spheres.push_back(sphere{vec3((radius - ball_radius * .7) * sin(i), current_height, (radius - ball_radius * .7) * cos(i)) + position, ball_radius, mat});
    }
}

void init_scene() {
    printf("init_scene\n");
    auto camera_pos = vec3(0,3,5) / 1.5;
    auto target_pos = vec3(-2, -1, 0);
    auto delta = target_pos - camera_pos;
    current_camera = camera(camera_pos, delta, v3_up, 30, bitmap_width * 1.0 / bitmap_height, 0.001, 3);

    auto diffuse_red =    new material(vec3(1,.2,.2), 0, 0, material_type::diffuse, false);
    auto diffuse_yellow = new material(vec3(.8,.8,.2), 0, 0, material_type::diffuse, false);
    auto diffuse_blue =   new material(vec3(.2,.2,.6), 0, 0, material_type::diffuse, false);
    auto diffuse_green =  new material(vec3(.2,.6,.2), 0, 0, material_type::diffuse, false);
    auto diffuse_white =  new material(vec3(.5,.5,.5), 0, 0, material_type::diffuse, false);
    auto diffuse_light =  new material(v3_one * 4, 0, 0, material_type::diffuse, true);

    auto chrome =         new material(v3_one * .9, 0, 0, material_type::metal, false);
    auto gold =           new material(vec3(.9,.7,.2), 0, 0, material_type::metal, false);
    auto aluminum =       new material(v3_one * .9, .01, 0, material_type::metal, false);
    auto glass =          new material(v3_one, 0, 1.5, material_type::dielectric, false);
    auto blue_gem =       new material(vec3(.1,.3,.7), 0, 2.4, material_type::dielectric, false);

    spheres.push_back(sphere{vec3(-3.5,1.5,0), 2, aluminum});
    spheres.push_back(sphere{vec3(0,-1003,-1), 1000, diffuse_blue});

    cones.push_back(cone{vec3(-1,1,0), -v3_up, 1.5, cos(M_PI / 8) * cos(M_PI / 8), diffuse_green});
    discs.push_back(disc{vec3(-1,-0.5,0), -v3_up, 1.5f * sin(M_PI / 8), diffuse_green});
    add_spiraling_spheres(vec3(-1,1,0), sin(M_PI / 6), 1.5, 15);

    cones.push_back(cone{vec3(-1,0,0), -v3_up, 1.5, cos(M_PI / 6) * cos(M_PI / 6), diffuse_green});
    discs.push_back(disc{vec3(-1,-1.5,0), -v3_up, 1.5f * sin(M_PI / 6), diffuse_blue});
    add_spiraling_spheres(vec3(-1,0,0), sin(M_PI / 4), 1.5, 20);

    cones.push_back(cone{vec3(-1,-1,0), -v3_up, 1.5, cos(M_PI / 5) * cos(M_PI / 5), diffuse_green});
    discs.push_back(disc{vec3(-1,-2.5,0), -v3_up, 1.5f * sin(M_PI / 5), diffuse_green});
    add_spiraling_spheres(vec3(-1,-1,0), sin(M_PI / 3), 1.5, 25);

    auto box_pos = vec3(-2.2, -2.7, 1);
    auto box_extends = vec3(1,1,.8) * .3;
    auto line1_extends = vec3(1.04,1.04, .2) * .3;
    auto line3_extends = vec3(.2,1.04, .84) * .3;
    auto line2_pos = box_pos + v3_up * box_extends.y * 1.1;
    auto line2_extends = .15 * v3_one;

    boxes.push_back(AABB{box_pos - box_extends, box_pos+box_extends, diffuse_red});
    boxes.push_back(AABB{box_pos - line1_extends, box_pos + line1_extends, diffuse_yellow});
    boxes.push_back(AABB{box_pos - line3_extends, box_pos + line3_extends, diffuse_yellow});
    boxes.push_back(AABB{line2_pos - line2_extends, line2_pos + line2_extends, diffuse_yellow});

    for (int i = 0; i < 12; i++) {
        auto pos = vec3(-1 + 10 * sin(2 * M_PI * i / 12), 1 +  0.1 * sin(i), 10 * cos(2 * M_PI * i / 12));
        spheres.push_back(sphere{pos, .5, diffuse_light});    
    }
}