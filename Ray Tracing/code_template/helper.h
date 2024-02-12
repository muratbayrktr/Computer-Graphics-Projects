#include <vector>
#include <iostream>
#include <cmath>
#include "parser.h"

void precompute_mesh_normals(const std::vector<parser::Mesh>& meshes, std::vector<std::vector<parser::Vec3f>>& mesh_normals, const std::vector<parser::Vec3f>& vertex_data);

void precompute_triangle_normals(const std::vector<parser::Triangle>& triangles, std::vector<parser::Vec3f>& triangle_normals, const std::vector<parser::Vec3f>& vertex_data);

parser::Ray cast_ray(const parser::Camera& camera, int x, int y, int width, int height);

parser::Vec3f compute_color(const parser::Scene& scene, parser::Ray& ray, const std::vector<std::vector<parser::Vec3f>>& mesh_normals, const std::vector<parser::Vec3f>& triangle_normals);

bool closest_hit(parser::Ray& ray, parser::hitRecord& hit_record, const parser::Scene& scene, const std::vector<std::vector<parser::Vec3f>>& mesh_normals, const std::vector<parser::Vec3f>& triangle_normals);

parser::Vec3f apply_shading(parser::Ray& ray, parser::hitRecord& hit_record, const parser::Scene& scene, const std::vector<std::vector<parser::Vec3f>>& mesh_normals, const std::vector<parser::Vec3f>& triangle_normals);

parser::Vec3f normalize(const parser::Vec3f& v);

parser::Vec3f cross_product(const parser::Vec3f& v1, const parser::Vec3f& v2);

float dot_product(const parser::Vec3f& v1, const parser::Vec3f& v2);

bool ray_face_intersection(const parser::Vec3f& ray_origin, const parser::Vec3f& ray_direction, 
                               const parser::Vec3f& A, const parser::Vec3f& B, const parser::Vec3f& C, 
                               float& t);

bool ray_sphere_intersection(const parser::Vec3f& ray_origin, const parser::Vec3f& ray_direction, 
                               const parser::Vec3f& center, float radius, parser::Vec3f& intersection_point, parser::Vec3f& intersection_normal, float& t_min);

bool ray_triangle_intersection(const parser::Vec3f& ray_origin, const parser::Vec3f& ray_direction, 
                           const std::vector<parser::Vec3f>& vertex_data, const parser::Triangle& triangle, 
                           const std::vector<parser::Vec3f>& normals, parser::Vec3f& intersection_point, 
                           parser::Vec3f& intersection_normal, float& t_min);

bool ray_mesh_intersection(const parser::Vec3f& ray_origin, const parser::Vec3f& ray_direction, 
                           const std::vector<parser::Vec3f>& vertex_data, const parser::Mesh& mesh, 
                           const std::vector<parser::Vec3f>& normals, parser::Vec3f& intersection_point, 
                           parser::Vec3f& intersection_normal, float& t_min);

bool is_in_shadow(const parser::Vec3f& point, const parser::Vec3f& intersection_normal, const parser::Vec3f& lightPosition, const parser::Scene& scene, const std::vector<std::vector<parser::Vec3f>>& mesh_normals, const std::vector<parser::Vec3f>& triangle_normals);

parser::Vec3f compute_diffuse(const parser::PointLight& light, const parser::Vec3f& intersection_point, 
                              const parser::Vec3f& intersection_normal, const parser::Vec3f& diffuse_reflectance, parser::Vec3f& diffuse_accumulator);

parser::Vec3f compute_specular(const parser::PointLight& light, const parser::Vec3f& intersection_point, 
                               const parser::Vec3f& intersection_normal, const parser::Vec3f& specular, 
                               const parser::Vec3f& ray_direction, float phong_exponent, parser::Vec3f& specular_accumulator);

parser::Ray reflect(parser::Ray& ray, const parser::Vec3f& intersection_normal, const parser::Vec3f& intersection_point, const float epsilon);                                                             

