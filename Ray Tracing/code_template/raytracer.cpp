#include <iostream>
#include "parser.h"
#include "ppm.h"
#include "helper.h"
#include <pthread.h>
typedef unsigned char RGB[3];

int width;
int height;
unsigned char* image;
parser::Scene scene;
pthread_t thread_1;
pthread_t thread_2;
pthread_t thread_3;
pthread_t thread_4;

typedef struct thread_args {
    int camera_id;
    int thread_id;
    std::vector<std::vector<parser::Vec3f>> mesh_normals;
    std::vector<parser::Vec3f> triangle_normals;
} thread_args;

void* compute_pixels(void * args) {
    thread_args *thread_args = (struct thread_args *) args;
    std::vector<std::vector<parser::Vec3f>> mesh_normals = thread_args->mesh_normals;
    std::vector<parser::Vec3f> triangle_normals = thread_args->triangle_normals;
    int i = thread_args->camera_id;

    int thread_id = thread_args->thread_id;

    for (int y = thread_id % 4; y < height; y += 4 ) {
        for (int x = 0; x < width; x++) {
            parser::Ray ray = cast_ray(scene.cameras[i], x, y, scene.cameras[i].image_width, scene.cameras[i].image_height);
            ray.depth = 0;
            parser::Vec3f image_color = compute_color(scene, ray, mesh_normals, triangle_normals);
            if (image_color.x > 255) {
                image_color.x = 255;
            }
            if (image_color.y > 255) {
                image_color.y = 255;
            }
            if (image_color.z > 255) {
                image_color.z = 255;
            }
            image[(y * width + x) * 3 + 0] = std::round(image_color.x);
            image[(y * width + x) * 3 + 1] = std::round(image_color.y);
            image[(y * width + x) * 3 + 2] = std::round(image_color.z);
        }
    }
    pthread_exit(NULL);
    return NULL;
}

int main(int argc, char* argv[])
{
    scene.loadFromXml(argv[1]);    
    std::vector<std::vector<parser::Vec3f>> mesh_normals;
    std::vector<parser::Vec3f> triangle_normals;
    precompute_mesh_normals(scene.meshes, mesh_normals, scene.vertex_data);
    precompute_triangle_normals(scene.triangles, triangle_normals, scene.vertex_data);
    for (int i = 0; i < scene.cameras.size(); i++) {
        width = scene.cameras[i].image_width;
        height = scene.cameras[i].image_height;
        image = new unsigned char [width * height * 3];
        thread_args args_1;
        args_1.camera_id = i;
        args_1.thread_id = 0;
        args_1.mesh_normals = mesh_normals;
        args_1.triangle_normals = triangle_normals;
        pthread_create(&thread_1, NULL, compute_pixels, (void *)&args_1);

        thread_args args_2;
        args_2.camera_id = i;
        args_2.thread_id = 1;
        args_2.mesh_normals = mesh_normals;
        args_2.triangle_normals = triangle_normals;
        pthread_create(&thread_2, NULL, compute_pixels, (void *)&args_2);

        thread_args args_3;
        args_3.camera_id = i;
        args_3.thread_id = 2;
        args_3.mesh_normals = mesh_normals;
        args_3.triangle_normals = triangle_normals;
        pthread_create(&thread_3, NULL, compute_pixels, (void *)&args_3);
        
        thread_args args_4;
        args_4.camera_id = i;
        args_4.thread_id = 3;
        args_4.mesh_normals = mesh_normals;
        args_4.triangle_normals = triangle_normals;
        pthread_create(&thread_4, NULL, compute_pixels, (void *)&args_4);


        pthread_join(thread_1, NULL);
        pthread_join(thread_2, NULL);
        pthread_join(thread_3, NULL);
        pthread_join(thread_4, NULL);

        write_ppm(scene.cameras[i].image_name.c_str(), image, width, height);
    }            
    delete [] image;
    return 0;
}
