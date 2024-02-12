#ifndef __HELPERS_H__
#define __HELPERS_H__
#define ABS(a) ((a) > 0 ? (a) : -1 * (a))
#define EPSILON 0.000000001
#include "Matrix4.h"
#include "Vec3.h"
#include "Vec4.h"
#include "Camera.h"
#include "Mesh.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Color.h"
#include <vector>
/*
 * Calculate cross product of vec3 a, vec3 b and return resulting vec3.
 */
Vec3 crossProductVec3(Vec3 a, Vec3 b);

/*
 * Calculate dot product of vec3 a, vec3 b and return resulting value.
 */
double dotProductVec3(Vec3 a, Vec3 b);

/*
 * Find length (|v|) of vec3 v.
 */
double magnitudeOfVec3(Vec3 v);

/*
 * Normalize the vec3 to make it unit vec3.
 */
Vec3 normalizeVec3(Vec3 v);

/*
 * Return -v (inverse of vec3 v)
 */
Vec3 inverseVec3(Vec3 v);

/*
 * Add vec3 a to vec3 b and return resulting vec3 (a+b).
 */
Vec3 addVec3(Vec3 a, Vec3 b);

/*
 * Subtract vec3 b from vec3 a and return resulting vec3 (a-b).
 */
Vec3 subtractVec3(Vec3 a, Vec3 b);

/*
 * Multiply each element of vec3 with scalar.
 */
Vec3 multiplyVec3WithScalar(Vec3 v, double c);

/*
 * Prints elements in a vec3. Can be used for debugging purposes.
 */
void printVec3(Vec3 v);

/*
 * Check whether vec3 a and vec3 b are equal.
 * In case of equality, returns 1.
 * Otherwise, returns 0.
 */
int areEqualVec3(Vec3 a, Vec3 b);
/*
 * Returns an identity matrix (values on the diagonal are 1, others are 0).
 */
Matrix4 getIdentityMatrix();

/*
 * Multiply matrices m1 (Matrix4) and m2 (Matrix4) and return the result matrix r (Matrix4).
 */
Matrix4 multiplyMatrixWithMatrix(Matrix4 m1, Matrix4 m2);

/*
 * Multiply matrix m (Matrix4) with vector v (vec4) and store the result in vector r (vec4).
 */
Vec4 multiplyMatrixWithVec4(Matrix4 m, Vec4 v);

Matrix4 getTranslationMatrix(double tx, double ty, double tz);

Matrix4 getScalingMatrix(double sx, double sy, double sz);

Matrix4 getRotationMatrix(Vec3 axis, double angle);

Matrix4 getCameraTransformation(Camera *camera);

Matrix4 getProjectionTransformation(Camera *camera);

Matrix4 getViewportTransformation(Camera *camera);

Matrix4 getModelingTransformation(Mesh *mesh, std::vector<Scaling *> scalings, std::vector<Rotation *> rotations, std::vector<Translation *> translations);

std::vector<Vec4 > getProjectedVertices(Triangle triangle, std::vector<Vec3 *> vertices, Matrix4 transformation);

std::vector<Vec4 > perspectiveDivision(std::vector<Vec4 > projectedVertices);

bool Culling(bool &cullingEnabled, std::vector<Vec4 > projectedVertices, Vec3 gaze);

bool CullingBefore(bool &cullingEnabled, std::vector<Vec4> projectedVertices, Vec3 gaze);

void triangleRasterization(std::vector<std::vector<Color>> & image, std::vector<std::vector<double> > &depth, std::vector<Color *> colorsOfVertices, std::vector<Vec4 > projectedVertices,  int nx, int ny);

bool visible(double den, double num, double &te, double &tl);

bool clippingLiangBarsky(Vec4& v0, Vec4& v1, Color& c0, Color& c1);

void lineRasterization(std::vector<std::vector<Color>> & image, std::vector<std::vector<double> > &depth, Vec4 v0, Vec4 v1, Color c0, Color c1);


#endif