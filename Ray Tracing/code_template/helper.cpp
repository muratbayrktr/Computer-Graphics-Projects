#include "helper.h"

// normalize
parser::Vec3f normalize(const parser::Vec3f& v) {
    parser::Vec3f result;
    float length = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    result.x = v.x / length;
    result.y = v.y / length;
    result.z = v.z / length;
    return result;
}

// Cross product 2 float vectors
parser::Vec3f cross_product(const parser::Vec3f& v1, const parser::Vec3f& v2) {
    parser::Vec3f result;
    result.x = v1.y * v2.z - v1.z * v2.y;
    result.y = v1.z * v2.x - v1.x * v2.z;
    result.z = v1.x * v2.y - v1.y * v2.x;
    return result;
}

// Dot product
float dot_product(const parser::Vec3f& v1, const parser::Vec3f& v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

// Precompute the normal vectors of each face of meshes
void precompute_mesh_normals(const std::vector<parser::Mesh>& meshes, 
                                std::vector<std::vector<parser::Vec3f>>& normals,
                                const std::vector<parser::Vec3f>& vertex_data) {
    normals.reserve(meshes.size());
    for (int i = 0; i < meshes.size(); i++) {
        normals.push_back(std::vector<parser::Vec3f>());
        for (int j = 0; j < meshes[i].faces.size(); j++) {
            parser::Vec3f v0 = vertex_data[meshes[i].faces[j].v0_id-1];
            parser::Vec3f v1 = vertex_data[meshes[i].faces[j].v1_id-1];
            parser::Vec3f v2 = vertex_data[meshes[i].faces[j].v2_id-1];
            // vectors are given counter-clockwise
            parser::Vec3f v1_v0, v2_v0;
            // minus(v1 , v0), minus(v2 , v0)
            v1_v0.x = v1.x - v0.x;
            v1_v0.y = v1.y - v0.y;
            v1_v0.z = v1.z - v0.z;
            v2_v0.x = v2.x - v0.x;
            v2_v0.y = v2.y - v0.y;
            v2_v0.z = v2.z - v0.z;

            parser::Vec3f normal = normalize(cross_product(v1_v0, v2_v0));
            normals[i].push_back(normal);
        }
    }
}
// Precompute the normal vectors of each triangle
void precompute_triangle_normals(const std::vector<parser::Triangle>& triangles, 
                                std::vector<parser::Vec3f>& normals,
                                const std::vector<parser::Vec3f>& vertex_data) {
    normals.reserve(triangles.size());
    for (int i = 0; i < triangles.size(); i++) {
        parser::Vec3f v0 = vertex_data[triangles[i].indices.v0_id-1];
        parser::Vec3f v1 = vertex_data[triangles[i].indices.v1_id-1];
        parser::Vec3f v2 = vertex_data[triangles[i].indices.v2_id-1];
        // vectors are given counter-clockwise
        parser::Vec3f v1_v0, v2_v0;
        // minus(v1 , v0), minus(v2 , v0)
        v1_v0.x = v1.x - v0.x;
        v1_v0.y = v1.y - v0.y;
        v1_v0.z = v1.z - v0.z;
        v2_v0.x = v2.x - v0.x;
        v2_v0.y = v2.y - v0.y;
        v2_v0.z = v2.z - v0.z;

        parser::Vec3f normal = normalize(cross_product(v1_v0, v2_v0));
        normals.push_back(normal);
    }
}


parser::Ray cast_ray(const parser::Camera& camera, int x, int y, int width, int height) {
    parser::Ray ray;
    ray.origin = camera.position;
    parser::Vec3f w = parser::Vec3f({-camera.gaze.x, -camera.gaze.y, -camera.gaze.z});
    w = normalize(w);
    parser::Vec3f v = camera.up;
    v = normalize(v);
    parser::Vec3f u = cross_product(v, w);
    u = normalize(u);
    float l = camera.near_plane.x;
    float r = camera.near_plane.y;
    float b = camera.near_plane.z;
    float t = camera.near_plane.w;
    float s_u = (x+0.5)*(r - l) / width;
    float s_v = (y+0.5)*(t - b) / height;

    // parser::Vec3f m = ray.origin + camera.gaze * camera.near_distance; 
    parser::Vec3f m, q, s, ray_direction;
    m.x = ray.origin.x + camera.gaze.x * camera.near_distance;
    m.y = ray.origin.y + camera.gaze.y * camera.near_distance;
    m.z = ray.origin.z + camera.gaze.z * camera.near_distance;

    // parser::Vec3f q = m + l*u + t*v; // top left corner
    q.x = m.x + l*u.x + t*v.x;
    q.y = m.y + l*u.y + t*v.y;
    q.z = m.z + l*u.z + t*v.z;

    // parser::Vec3f s = q + u*s_u - v*s_v; // corresponding point on the image plane
    s.x = q.x + u.x*s_u - v.x*s_v;
    s.y = q.y + u.y*s_u - v.y*s_v;
    s.z = q.z + u.z*s_u - v.z*s_v;

    // parser::Vec3f ray_direction = s - ray.origin; // ray direction
    ray_direction.x = s.x - ray.origin.x;
    ray_direction.y = s.y - ray.origin.y;
    ray_direction.z = s.z - ray.origin.z;

    // normalize the ray direction
    ray_direction = normalize(ray_direction); 
    ray.direction = ray_direction;
    return ray;
}

bool ray_face_intersection(const parser::Vec3f& ray_origin, const parser::Vec3f& ray_direction, 
                               const parser::Vec3f& A, const parser::Vec3f& B, const parser::Vec3f& C, 
                               float& t) {
    parser::Vec3f AB = parser::Vec3f({B.x - A.x, B.y - A.y, B.z - A.z});
    parser::Vec3f AC = parser::Vec3f({C.x - A.x, C.y - A.y, C.z - A.z});

    parser::Vec3f h = cross_product(ray_direction,AC);
    float a = dot_product(AB,h);
    float my_epsilon = 1e-20;

    if (a > -my_epsilon && a < my_epsilon)
        return false;  // Ray is parallel to the triangle

    float f = 1.0 / a;
    parser::Vec3f s = parser::Vec3f({ray_origin.x - A.x, ray_origin.y - A.y, ray_origin.z - A.z});
    float u = f * dot_product(s,h);

    if (u < 0.0 || u > 1.0)
        return false; 

    parser::Vec3f q = cross_product(s,AB);
    float v = f * dot_product(ray_direction,q);

    if (v < 0.0 || u + v > 1.0)
        return false;

    t = f * dot_product(AC,q); 

    if (t > my_epsilon)  // ray intersection
        return true; 

    return false;  // line intersection but not ray intersection
}

bool ray_sphere_intersection(const parser::Vec3f& ray_origin, const parser::Vec3f& ray_direction, 
                               const parser::Vec3f& center, float radius, parser::Vec3f& intersection_point, parser::Vec3f& intersection_normal, float& t_min) {
    
    // Vector from the sphere's center to the ray's origin
    parser::Vec3f oc = {ray_origin.x - center.x, ray_origin.y - center.y, ray_origin.z - center.z};
    // Coefficients for the quadratic formula
    float a = dot_product(ray_direction, ray_direction);
    float b = 2.0f * dot_product(oc, ray_direction);
    float c_val = dot_product(oc, oc) - radius * radius;

    // Discriminant of the quadratic
    float discriminant = b * b - 4 * a * c_val;

    // If the discriminant is negative, no intersection
    if (discriminant < 0) {
        return false;
    } else {
        // Closest intersection point
        float t1 = 0.5 * (-b - sqrt(discriminant)) / a;
        float t2 = 0.5 * (-b + sqrt(discriminant)) / a;
        if((t1 <= t_min) && (t1 <= t2) && (t1 > 0))
        {
            t_min = t1;
            intersection_point = {
            ray_origin.x + t1 * ray_direction.x,
            ray_origin.y + t1 * ray_direction.y,
            ray_origin.z + t1 * ray_direction.z
            };
            intersection_normal.x = (intersection_point.x - center.x) / radius;
            intersection_normal.y = (intersection_point.y - center.y) / radius;
            intersection_normal.z = (intersection_point.z - center.z) / radius;
            return true;
        }
        else if((t2 < t_min) && (t2 < t1) && t2 > 0)
        {
            t_min = t2;
            intersection_point = {
            ray_origin.x + t2 * ray_direction.x,
            ray_origin.y + t2 * ray_direction.y,
            ray_origin.z + t2 * ray_direction.z
            };
            intersection_normal.x = (intersection_point.x - center.x) / radius;
            intersection_normal.y = (intersection_point.y - center.y) / radius;
            intersection_normal.z = (intersection_point.z - center.z) / radius;
            return true;
        }
        return false;
    }
}

bool ray_triangle_intersection(const parser::Vec3f& ray_origin, const parser::Vec3f& ray_direction, 
                           const std::vector<parser::Vec3f>& vertex_data, const parser::Triangle& triangle, 
                           const std::vector<parser::Vec3f>& normals, parser::Vec3f& intersection_point, 
                           parser::Vec3f& intersection_normal, float& t_min) {
    const parser::Vec3f& A = vertex_data[triangle.indices.v0_id-1];
    const parser::Vec3f& B = vertex_data[triangle.indices.v1_id-1];
    const parser::Vec3f& C = vertex_data[triangle.indices.v2_id-1];
    if(dot_product(normals[triangle.material_id-1], ray_direction) > 0)
    {
        return false;
    }
    float t;
    if (ray_face_intersection(ray_origin, ray_direction, A, B, C, t)) {
        if (t < t_min) {
            t_min = t;
            // add(ray_origin, times(t, ray_direction));
            intersection_point = parser::Vec3f({ray_origin.x + t * ray_direction.x, 
                        ray_origin.y + t * ray_direction.y, ray_origin.z + t * ray_direction.z});
            intersection_normal = normals[triangle.material_id-1];
            // intersection_normal = cross_product(B - A,C - A);  // assuming counter-clockwise winding
            // intersection_normal = normalize(intersection_normal);  // normalize the normal
        }
        return true;
    }
    return false;
}

bool ray_mesh_intersection(const parser::Vec3f& ray_origin, const parser::Vec3f& ray_direction, 
                           const std::vector<parser::Vec3f>& vertex_data, const parser::Mesh& mesh, 
                           const std::vector<parser::Vec3f>& normals, parser::Vec3f& intersection_point, 
                           parser::Vec3f& intersection_normal, float& t_min) {
    bool intersected = false;
    
    for (int i = 0; i < mesh.faces.size(); i++) {
        const auto& face = mesh.faces[i];
        const parser::Vec3f& A = vertex_data[face.v0_id-1];
        const parser::Vec3f& B = vertex_data[face.v1_id-1];
        const parser::Vec3f& C = vertex_data[face.v2_id-1];

        float t;
        if (ray_face_intersection(ray_origin, ray_direction, A, B, C, t)) {
            if (t < t_min) {
                t_min = t;
                intersected = true;
                // add(ray_origin, times(t, ray_direction));
                intersection_point = parser::Vec3f({ray_origin.x + t * ray_direction.x, 
                        ray_origin.y + t * ray_direction.y, ray_origin.z + t * ray_direction.z});
                intersection_normal = normals[i];
                // intersection_normal = cross_product(B - A,C - A);  // assuming counter-clockwise winding
                // intersection_normal = normalize(intersection_normal);  // normalize the normal
            }
        }
    }
    return intersected;
}

bool closest_hit(parser::Ray& ray, parser::hitRecord& hit_record, const parser::Scene& scene, const std::vector<std::vector<parser::Vec3f>>& mesh_normals, const std::vector<parser::Vec3f>& triangle_normals){
    float t_min = 3.40282e+038;
    bool hit = false;
    hit_record.hit = false;
    // for each object in the scene check if the ray intersects with it and get t min
    for (int i=0; i < scene.triangles.size();i++) {
        parser::Vec3f intersection_point, intersection_normal;
        hit = ray_triangle_intersection(ray.origin, ray.direction, scene.vertex_data, scene.triangles[i], triangle_normals, intersection_point, intersection_normal, t_min);
        if (hit) {
            hit_record.intersection_point = intersection_point;
            hit_record.intersection_normal = intersection_normal;
            hit_record.material = scene.materials[scene.triangles[i].material_id-1];
            hit_record.hit = true;

        }
    }
    for (int i=0; i < scene.meshes.size();i++) {
        parser::Vec3f intersection_point, intersection_normal;
        hit = ray_mesh_intersection(ray.origin, ray.direction, scene.vertex_data, scene.meshes[i], mesh_normals[i], intersection_point, intersection_normal, t_min);
        if (hit) {
            hit_record.intersection_point = intersection_point;
            hit_record.intersection_normal = intersection_normal;
            hit_record.material = scene.materials[scene.meshes[i].material_id-1];
            hit_record.hit = true;

        }
    }
    for (int i=0; i < scene.spheres.size();i++) {
        parser::Vec3f intersection_point, intersection_normal;
        hit = ray_sphere_intersection(ray.origin, ray.direction, scene.vertex_data[scene.spheres[i].center_vertex_id-1], scene.spheres[i].radius, intersection_point, intersection_normal, t_min);
        if (hit) {
            hit_record.intersection_point = intersection_point;
            hit_record.intersection_normal = intersection_normal;
            hit_record.material = scene.materials[scene.spheres[i].material_id-1];
            hit_record.hit = true;
        }
    }
    return hit_record.hit;
}


bool is_in_shadow(const parser::Vec3f& point, const parser::Vec3f& intersection_normal, const parser::Vec3f& lightPosition, const parser::Scene& scene, const std::vector<std::vector<parser::Vec3f>>& mesh_normals, const std::vector<parser::Vec3f>& triangle_normals) {
    parser::Vec3f toLight = parser::Vec3f({lightPosition.x - point.x, lightPosition.y - point.y, lightPosition.z - point.z});
    float distanceToLight = sqrt(dot_product(toLight, toLight));
    toLight = normalize(toLight);

    float tNearShadow = 3.40282e+038;

    // Check meshes for shadow intersection
    for (int m = 0; m < scene.meshes.size(); m++) {
        parser::Vec3f shadowIntersection = {0, 0, 0};
        parser::Vec3f shadowNormal = mesh_normals[m][0];

        if (ray_mesh_intersection(parser::Vec3f({point.x + scene.shadow_ray_epsilon * intersection_normal.x, 
                                                 point.y + scene.shadow_ray_epsilon * intersection_normal.y, 
                                                 point.z + scene.shadow_ray_epsilon * intersection_normal.z}),
                                  toLight, 
                                  scene.vertex_data, 
                                  scene.meshes[m], 
                                  mesh_normals[m], 
                                  shadowIntersection, 
                                  shadowNormal, 
                                  tNearShadow) && tNearShadow < distanceToLight) {
            
            return true;
        }
    }

    // Check triangles for shadow intersection
    for (int m = 0; m < scene.triangles.size(); m++) {
        parser::Vec3f shadowIntersection = {0, 0, 0};
        parser::Vec3f shadowNormal = triangle_normals[m];

        if (ray_triangle_intersection(parser::Vec3f({point.x + scene.shadow_ray_epsilon * intersection_normal.x, 
                                                 point.y + scene.shadow_ray_epsilon * intersection_normal.y, 
                                                 point.z + scene.shadow_ray_epsilon * intersection_normal.z}),
                                      toLight, 
                                      scene.vertex_data, 
                                      scene.triangles[m], 
                                      triangle_normals, 
                                      shadowIntersection, 
                                      shadowNormal, 
                                      tNearShadow) && tNearShadow < distanceToLight) {
            
            return true;
        }
    }

    // Check spheres for shadow intersection
    for (int m = 0; m < scene.spheres.size(); m++) {
        parser::Vec3f shadowIntersection = {0, 0, 0};
        parser::Vec3f shadowNormal = {0, 0, 0};
        if (ray_sphere_intersection(parser::Vec3f({point.x + scene.shadow_ray_epsilon * intersection_normal.x, 
                                                 point.y + scene.shadow_ray_epsilon * intersection_normal.y, 
                                                 point.z + scene.shadow_ray_epsilon * intersection_normal.z}),
                                    toLight, 
                                    scene.vertex_data[scene.spheres[m].center_vertex_id-1],
                                    scene.spheres[m].radius,
                                    shadowIntersection, shadowNormal, 
                                    tNearShadow) && tNearShadow < distanceToLight) {
            
            return true;
        }
    }
    return false;
}

parser::Vec3f compute_diffuse(const parser::PointLight& light, const parser::Vec3f& intersection_point, 
                              const parser::Vec3f& intersection_normal, const parser::Vec3f& diffuse_reflectance, parser::Vec3f& diffuse_accumulator) {
    
    // Compute the vector from the intersection point to the light source.
    parser::Vec3f light_vector = parser::Vec3f({light.position.x - intersection_point.x, 
                                                light.position.y - intersection_point.y, 
                                                light.position.z - intersection_point.z});

    // Calculate the magnitude of this vector which gives the distance 'r' from the light source.
    float r = sqrt(light_vector.x * light_vector.x + light_vector.y * light_vector.y + light_vector.z * light_vector.z);

    // Normalize the light vector.
    light_vector.x /= r;
    light_vector.y /= r;
    light_vector.z /= r;

    // Compute the dot product between this vector and the intersection normal.
    float dot_product = std::max(0.0f, light_vector.x * intersection_normal.x + light_vector.y * intersection_normal.y + light_vector.z * intersection_normal.z);

    // Calculate the adjusted intensity based on the inverse square law.
    parser::Vec3f adjusted_intensity;
    adjusted_intensity.x = light.intensity.x / (r * r);
    adjusted_intensity.y = light.intensity.y / (r * r);
    adjusted_intensity.z = light.intensity.z / (r * r);

    // Update the accumulator using the adjusted intensity.
    
    diffuse_accumulator.x += dot_product * adjusted_intensity.x * diffuse_reflectance.x;
    diffuse_accumulator.y += dot_product * adjusted_intensity.y * diffuse_reflectance.y;
    diffuse_accumulator.z += dot_product * adjusted_intensity.z * diffuse_reflectance.z;

    return diffuse_accumulator; 
}

parser::Vec3f compute_specular(const parser::PointLight& light, const parser::Vec3f& intersection_point, 
                               const parser::Vec3f& intersection_normal, const parser::Vec3f& specular, 
                               const parser::Vec3f& ray_direction, float phong_exponent, parser::Vec3f& specular_accumulator) {

    parser::Vec3f l = parser::Vec3f({light.position.x - intersection_point.x, 
                                     light.position.y - intersection_point.y, 
                                     light.position.z - intersection_point.z});
    parser::Vec3f l_norm = normalize(l);
    if(dot_product(l_norm, intersection_normal) >= 90)
    {
        parser::Vec3f specular_color = parser::Vec3f({0, 0, 0});
        return specular_color;
    }
    parser::Vec3f v = parser::Vec3f({-ray_direction.x, -ray_direction.y, -ray_direction.z});
    parser::Vec3f h = parser::Vec3f({l_norm.x + v.x, l_norm.y + v.y, l_norm.z + v.z});
    h = normalize(h);

    float cos_alpha = std::max(0.0F, dot_product(h, intersection_normal));
    // Calculate the magnitude of this vector which gives the distance 'r' from the light source.
    float r = sqrt(l.x * l.x + l.y * l.y + l.z * l.z);

    // Calculate the adjusted intensity based on the inverse square law.
    parser::Vec3f adjusted_intensity;
    adjusted_intensity.x = light.intensity.x / (r * r);
    adjusted_intensity.y = light.intensity.y / (r * r);
    adjusted_intensity.z = light.intensity.z / (r * r);


    parser::Vec3f specular_color; // specular * adjusted_intensity;
    specular_color.x = specular.x * adjusted_intensity.x;
    specular_color.y = specular.y * adjusted_intensity.y;
    specular_color.z = specular.z * adjusted_intensity.z;
    float coeff = pow(cos_alpha, phong_exponent);
    // Check for floating point errors
    
    specular_color = parser::Vec3f({specular_color.x * coeff, specular_color.y * coeff, specular_color.z * coeff});

    specular_accumulator.x += specular_color.x;
    specular_accumulator.y += specular_color.y;
    specular_accumulator.z += specular_color.z;

    return specular_color; // Not important we use accumulator
}
parser::Ray reflect(parser::Ray& ray, const parser::Vec3f& intersection_normal, const parser::Vec3f& intersection_point, const float epsilon) {
    parser::Ray reflected_ray;
    reflected_ray.origin = ray.origin;
    reflected_ray.direction = ray.direction;
    reflected_ray.depth = ray.depth;
    // Reflect the ray by changing the values of ray.direction and ray.origin
    reflected_ray.depth++;
    // wr = -wo + 2ncosΘ = -wo + 2n(n.wo)
    reflected_ray.direction.x = ray.direction.x - 2 * dot_product(ray.direction, intersection_normal) * intersection_normal.x;
    reflected_ray.direction.y = ray.direction.y - 2 * dot_product(ray.direction, intersection_normal) * intersection_normal.y;
    reflected_ray.direction.z = ray.direction.z - 2 * dot_product(ray.direction, intersection_normal) * intersection_normal.z;
    reflected_ray.origin.x = intersection_point.x - epsilon * ray.direction.x;
    reflected_ray.origin.y = intersection_point.y - epsilon * ray.direction.y;
    reflected_ray.origin.z = intersection_point.z - epsilon * ray.direction.z;

    return reflected_ray;
}

parser::Vec3f apply_shading(parser::Ray& ray, parser::hitRecord& hit_record, const parser::Scene& scene, const std::vector<std::vector<parser::Vec3f>>& mesh_normals, const std::vector<parser::Vec3f>& triangle_normals){
    parser::Material m = hit_record.material;
    parser::Vec3f mirror_reflectance = m.mirror;
    parser::Vec3f color = parser::Vec3f({0, 0, 0});
    // Ambient Part
    color.x += m.ambient.x * scene.ambient_light.x;
    color.y += m.ambient.y * scene.ambient_light.y;
    color.z += m.ambient.z * scene.ambient_light.z;
    // Mirror Part
    if (m.is_mirror) {
        // static parser::Ray reflected_ray;
        // reflected_ray.origin = ray.origin;
        // reflected_ray.direction = ray.direction;
        // reflected_ray.depth = ray.depth;
        parser::Ray reflected_ray = reflect(ray, hit_record.intersection_normal, hit_record.intersection_point, scene.shadow_ray_epsilon);
        // std::cout << "depth: " << ray.depth << std::endl;
        parser::Vec3f mirror_color = compute_color(scene, reflected_ray, mesh_normals, triangle_normals);
        color.x += mirror_reflectance.x * mirror_color.x;
        color.y += mirror_reflectance.y * mirror_color.y;
        color.z += mirror_reflectance.z * mirror_color.z;    
        // std::cout << "mirror: " << mirror_color.x << " " << mirror_color.y << " " << mirror_color.z << std::endl;
    }
    for(int i=0; i<scene.point_lights.size(); i++) {
        parser::PointLight light = scene.point_lights[i];
        if (!is_in_shadow(hit_record.intersection_point, hit_record.intersection_normal, light.position, scene, mesh_normals, triangle_normals)) {
            parser::Vec3f diffuseTerm, specularTerm;
            diffuseTerm.x = 0; diffuseTerm.y = 0; diffuseTerm.z = 0;
            specularTerm.x = 0; specularTerm.y = 0; specularTerm.z = 0;
            compute_diffuse(light, hit_record.intersection_point, hit_record.intersection_normal, m.diffuse, diffuseTerm);
            compute_specular(light, hit_record.intersection_point, hit_record.intersection_normal, m.specular, ray.direction, m.phong_exponent, specularTerm);
            color.x += diffuseTerm.x + specularTerm.x;
            color.y += diffuseTerm.y + specularTerm.y;
            color.z += diffuseTerm.z + specularTerm.z;
        }
    }
    return color;
}

parser::Vec3f compute_color(const parser::Scene& scene, parser::Ray& ray, const std::vector<std::vector<parser::Vec3f>>& mesh_normals, const std::vector<parser::Vec3f>& triangle_normals){
    if(ray.depth > scene.max_recursion_depth)
    {
        // std::cout << "Max recursion depth reached" << std::endl;
        return {0, 0, 0};
    }
    parser::hitRecord hit_record;
    if(closest_hit(ray, hit_record, scene, mesh_normals, triangle_normals))
    {
        return apply_shading(ray, hit_record, scene, mesh_normals, triangle_normals);
    }
    else if(ray.depth == 0)
    {
        parser::Vec3f bg;
        bg.x = scene.background_color.x;
        bg.y = scene.background_color.y;
        bg.z = scene.background_color.z;
        return bg;
    }
    return {0, 0, 0};
}
