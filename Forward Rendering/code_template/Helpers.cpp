#include <iostream>
#include <cmath>
#include "Helpers.h"


/*
 * Calculate cross product of vec3 a, vec3 b and return resulting vec3.
 */
Vec3 crossProductVec3(Vec3 a, Vec3 b)
{
    return Vec3(a.y * b.z - b.y * a.z, b.x * a.z - a.x * b.z, a.x * b.y - b.x * a.y);
}

/*
 * Calculate dot product of vec3 a, vec3 b and return resulting value.
 */
double dotProductVec3(Vec3 a, Vec3 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/*
 * Find length (|v|) of vec3 v.
 */
double magnitudeOfVec3(Vec3 v)
{
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

/*
 * Normalize the vec3 to make it unit vec3.
 */
Vec3 normalizeVec3(Vec3 v)
{
    double d = magnitudeOfVec3(v);
    return Vec3(v.x / d, v.y / d, v.z / d);
}

/*
 * Return -v (inverse of vec3 v)
 */
Vec3 inverseVec3(Vec3 v)
{
    return Vec3(-v.x, -v.y, -v.z);
}

/*
 * Add vec3 a to vec3 b and return resulting vec3 (a+b).
 */
Vec3 addVec3(Vec3 a, Vec3 b)
{
    return Vec3(a.x + b.x, a.y + b.y, a.z + b.z);
}

/*
 * Subtract vec3 b from vec3 a and return resulting vec3 (a-b).
 */
Vec3 subtractVec3(Vec3 a, Vec3 b)
{
    return Vec3(a.x - b.x, a.y - b.y, a.z - b.z);
}

/*
 * Multiply each element of vec3 with scalar.
 */
Vec3 multiplyVec3WithScalar(Vec3 v, double c)
{
    return Vec3(v.x * c, v.y * c, v.z * c);
}

/*
 * Prints elements in a vec3. Can be used for debugging purposes.
 */
void printVec3(Vec3 v)
{
    std::cout << "(" << v.x << "," << v.y << "," << v.z << ")" << std::endl;
}

/*
 * Check whether vec3 a and vec3 b are equal.
 * In case of equality, returns 1.
 * Otherwise, returns 0.
 */
int areEqualVec3(Vec3 a, Vec3 b)
{

    /* if x difference, y difference and z difference is smaller than threshold, then they are equal */
    if ((ABS((a.x - b.x)) < EPSILON) && (ABS((a.y - b.y)) < EPSILON) && (ABS((a.z - b.z)) < EPSILON))
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

/*
 * Returns an identity matrix (values on the diagonal are 1, others are 0).
 */
Matrix4 getIdentityMatrix()
{
    Matrix4 result;

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            if (i == j)
            {
                result.values[i][j] = 1.0;
            }
            else
            {
                result.values[i][j] = 0.0;
            }
        }
    }

    return result;
}

/*
 * Multiply matrices m1 (Matrix4) and m2 (Matrix4) and return the result matrix r (Matrix4).
 */
Matrix4 multiplyMatrixWithMatrix(Matrix4 m1, Matrix4 m2)
{
    Matrix4 result;
    double total;

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            total = 0;
            for (int k = 0; k < 4; k++)
            {
                total += m1.values[i][k] * m2.values[k][j];
            }

            result.values[i][j] = total;
        }
    }

    return result;
}

/*
 * Multiply matrix m (Matrix4) with vector v (vec4) and store the result in vector r (vec4).
 */
Vec4 multiplyMatrixWithVec4(Matrix4 m, Vec4 v)
{
    double values[4];
    double total;

    for (int i = 0; i < 4; i++)
    {
        total = 0;
        for (int j = 0; j < 4; j++)
        {
            total += m.values[i][j] * v.getNthComponent(j);
        }
        values[i] = total;
    }

    return Vec4(values[0], values[1], values[2], values[3], v.colorId);
}
// DOUBLE CHECKED, SHOULD BE CORRECT
Matrix4 getTranslationMatrix(double tx, double ty, double tz)
{
    Matrix4 result = getIdentityMatrix();
    result.values[0][3] = tx;
    result.values[1][3] = ty;
    result.values[2][3] = tz;
    return result;
}
// DOUBLE CHECKED, SHOULD BE CORRECT
Matrix4 getScalingMatrix(double sx, double sy, double sz)
{
    Matrix4 result = getIdentityMatrix();
    result.values[0][0] = sx;
    result.values[1][1] = sy;
    result.values[2][2] = sz;
    return result;
}
// DOUBLE CHECKED, SHOULD BE CORRECT
Matrix4 getRotationMatrix(Vec3 axis, double angle)
{
    // Input file given as ux, uy, uz. Directly implement Rotation About Arbitrary Axis using Alternative Method from Slides;
    Matrix4 result = getIdentityMatrix();
    /*  To find v, set the smallest component of u (in an absolute sense) to zero and
        swap the other two while
        negating one */
    
    Vec3 u = normalizeVec3(axis);
    Vec3 v;
    if (ABS(u.x) <= ABS(u.y) && ABS(u.x) <= ABS(u.z))
    {
        v = Vec3(0.0, -u.z, u.y);
    }
    else if (ABS(u.y) <= ABS(u.x) && ABS(u.y) <= ABS(u.z))
    {
        v = Vec3(-u.z, 0.0, u.x);
    }
    else
    {
        v = Vec3(-u.y, u.x, 0.0);
    }
    // w = u x v
    Vec3 w = crossProductVec3(u, v);
    // Normalize v and w
    v = normalizeVec3(v);
    w = normalizeVec3(w);
    // We assumed  that the origin of uvw is the same as the origin of xyz, 
    // lets translate origin of uvw to origin of xyz to handle this case
    Matrix4 T = getTranslationMatrix(-u.x, -u.y, -u.z);
    Matrix4 TInverse = getTranslationMatrix(u.x, u.y, u.z);
    // Now rotate uvw such that it aligns with xyz: call this transform M
    Matrix4 M = getIdentityMatrix();
    M.values[0][0] = u.x;
    M.values[0][1] = u.y;
    M.values[0][2] = u.z;
    M.values[1][0] = v.x;
    M.values[1][1] = v.y;
    M.values[1][2] = v.z;
    M.values[2][0] = w.x;
    M.values[2][1] = w.y;
    M.values[2][2] = w.z;
    // Now rotate about x axis by angle
    Matrix4 R = getIdentityMatrix();
    R.values[1][1] = cos(angle * M_PI / 180.0);
    R.values[1][2] = -sin(angle * M_PI / 180.0);
    R.values[2][1] = sin(angle * M_PI / 180.0);
    R.values[2][2] = cos(angle * M_PI / 180.0);
    // Undo the initial rotation: call this M^-1
    Matrix4 MInverse = getIdentityMatrix();
    MInverse.values[0][0] = u.x;
    MInverse.values[1][0] = u.y;
    MInverse.values[2][0] = u.z;
    MInverse.values[0][1] = v.x;
    MInverse.values[1][1] = v.y;
    MInverse.values[2][1] = v.z;
    MInverse.values[0][2] = w.x;
    MInverse.values[1][2] = w.y;
    MInverse.values[2][2] = w.z;
    // The final rotation transform is: T^-1 * M^-1 * R * M * T
    result = multiplyMatrixWithMatrix(M, T);
    result = multiplyMatrixWithMatrix(R, result);
    result = multiplyMatrixWithMatrix(MInverse, result);
    result = multiplyMatrixWithMatrix(TInverse, result);
    return result;
}
// DOUBLE CHECKED, SHOULD BE CORRECT
Matrix4 getCameraTransformation(Camera *camera) {
    Matrix4 result = getIdentityMatrix();
    Matrix4 T = getTranslationMatrix(-(camera->position.x), -(camera->position.y), -(camera->position.z));
    Matrix4 M = getIdentityMatrix();  // Rotation matrix
    
    // Camera's u, v, and w vectors should form the columns of M
    Vec3 u = normalizeVec3(camera->u); // Assuming these are already normalized
    Vec3 v = normalizeVec3(camera->v);
    Vec3 w = normalizeVec3(camera->w);
    
    M.values[0][0] = u.x; M.values[0][1] = u.y; M.values[0][2] = u.z;
    M.values[1][0] = v.x; M.values[1][1] = v.y; M.values[1][2] = v.z;
    M.values[2][0] = w.x; M.values[2][1] = w.y; M.values[2][2] = w.z;

    // First apply rotation (M), then translation (T)
    result = multiplyMatrixWithMatrix(M, T);
    return result;
}

// DOUBLE CHECKED, SHOULD BE CORRECT
Matrix4 getProjectionTransformation(Camera *camera)
{
    Matrix4 result;
    if(camera->projectionType == ORTOGRAPHIC_PROJECTION) // do orthographic projection
    {
        result = getIdentityMatrix();
        result.values[0][0] = 2.0 / (camera->right - camera->left);
        result.values[1][1] = 2.0 / (camera->top - camera->bottom);
        result.values[2][2] = -(2.0 / (camera->far - camera->near));
        result.values[0][3] = -((camera->right + camera->left) / (camera->right - camera->left));
        result.values[1][3] = -((camera->top + camera->bottom) / (camera->top - camera->bottom));
        result.values[2][3] = -((camera->far + camera->near) / (camera->far - camera->near));
    }
    else if (camera->projectionType == PERSPECTIVE_PROJECTION)// do perspective projection
    {
        result = getIdentityMatrix();
        result.values[0][0] = (2.0 * camera->near) / (camera->right - camera->left);
        result.values[1][1] = (2.0 * camera->near) / (camera->top - camera->bottom);
        result.values[2][2] = -((camera->far + camera->near) / (camera->far - camera->near));
        result.values[2][3] = -((2.0 * camera->near * camera->far) / (camera->far - camera->near));
        result.values[3][2] = -1.0;
        result.values[3][3] = 0.0;
        result.values[0][2] = (camera->right + camera->left) / (camera->right - camera->left);
        result.values[1][2] = (camera->top + camera->bottom) / (camera->top - camera->bottom);
    }
    return result;
}
//  DOUBLE CHECKED, SHOULD BE CORRECT
Matrix4 getViewportTransformation(Camera *camera)
{
    Matrix4 result = getIdentityMatrix();
    result.values[0][0] = camera->horRes / 2.0;
    result.values[1][1] = camera->verRes / 2.0;
    result.values[2][2] = 0.5;
    result.values[0][3] = (camera->horRes-1) / 2.0;
    result.values[1][3] = (camera->verRes-1) / 2.0;
    result.values[2][3] = 0.5;
    return result;
}
// DOUBLE CHECKED, SHOULD BE CORRECT
Matrix4 getModelingTransformation(Mesh *mesh, std::vector<Scaling *> scalings, std::vector<Rotation *> rotations, std::vector<Translation *> translations)
{
    Matrix4 result = getIdentityMatrix();
    // Do the transformations in the order
    for(int i = 0; i < mesh->numberOfTransformations; ++i)
    {
        if(mesh->transformationTypes[i] == 't')
        {
            Translation* t = translations[mesh->transformationIds[i]-1];
            result = multiplyMatrixWithMatrix(getTranslationMatrix(t->tx, t->ty, t->tz), result);
        }
        else if(mesh->transformationTypes[i] == 's')
        {
            Scaling* s = scalings[mesh->transformationIds[i]-1];
            result = multiplyMatrixWithMatrix(getScalingMatrix(s->sx, s->sy, s->sz), result);
        }
        else if(mesh->transformationTypes[i] == 'r')
        {
            Rotation* r = rotations[mesh->transformationIds[i]-1];
            result = multiplyMatrixWithMatrix(getRotationMatrix(Vec3(r->ux, r->uy, r->uz, -1.0), r->angle), result);
        }
    }
    return result;
}
// DOUBLE CHECKED, SHOULD BE CORRECT
std::vector<Vec4 > getProjectedVertices(Triangle triangle, std::vector<Vec3 *> vertices, Matrix4 transformation)
{
    std::vector<Vec4 > result;
    for(int i = 0; i < 3; i++)
    {
        Vec4 v = Vec4(vertices[triangle.vertexIds[i]-1]->x, vertices[triangle.vertexIds[i]-1]->y, vertices[triangle.vertexIds[i]-1]->z, 1.0, vertices[triangle.vertexIds[i]-1]->colorId);
        v = multiplyMatrixWithVec4(transformation, v);
        result.push_back(v);
    }
    return result;
}
// DOUBLE CHECKED, SHOULD BE CORRECT
bool Culling(bool &cullingEnabled, std::vector<Vec4> projectedVertices, Vec3 gaze)
{
    if(cullingEnabled)
    {
        Vec3 v1 = Vec3(projectedVertices[0].x, projectedVertices[0].y, projectedVertices[0].z);
        Vec3 v2 = Vec3(projectedVertices[1].x, projectedVertices[1].y, projectedVertices[1].z);
        Vec3 v3 = Vec3(projectedVertices[2].x, projectedVertices[2].y, projectedVertices[2].z);
        Vec3 v1v2 = subtractVec3(v1, v2);
        Vec3 v1v3 = subtractVec3(v3, v1);
        //Vec3 normal = {v1v2.y*v1v3.z - v1v2.z*v1v3.y, v1v2.z*v1v3.x - v1v2.x*v1v3.z, v1v2.x*v1v3.y - v1v2.y*v1v3.x};
        Vec3 normal = normalizeVec3(crossProductVec3(v1v2, v1v3));
        double result = dotProductVec3(normal, v1);
        if(result > 0.0)
        {
            return true; // backFace from slides ??
        }
        else
        {
            return false; // frontFace from slides ??
        }
    }
    else
    {
        return false;
    }
}
// DOUBLE CHECKED, SHOULD BE CORRECT
std::vector<Vec4 > perspectiveDivision(std::vector<Vec4 > projectedVertices)
{
    std::vector<Vec4 > result;
    for(int i = 0; i < projectedVertices.size(); i++)
    {
        Vec4 v = Vec4(projectedVertices[i].x / projectedVertices[i].t, projectedVertices[i].y / projectedVertices[i].t, projectedVertices[i].z / projectedVertices[i].t, 1.0, projectedVertices[i].colorId);
        result.push_back(v);
    }
    return result;
}
// DOUBLE CHECKED, SHOULD BE CORRECT
void triangleRasterization(std::vector<std::vector<Color>> &image, std::vector<std::vector<double> > &depth, std::vector<Color *> colorsOfVertices, std::vector<Vec4 > projectedVertices,  int nx, int ny)
{
    // Get vertices and its colors
    Vec4 v0 = projectedVertices[0];
    Vec4 v1 = projectedVertices[1];
    Vec4 v2 = projectedVertices[2];
    Color* c0 = colorsOfVertices[v0.colorId-1];
    Color* c1 = colorsOfVertices[v1.colorId-1];
    Color* c2 = colorsOfVertices[v2.colorId-1];
    
    // Get bounding box
    int xMin = std::min(std::min(v0.x, v1.x), v2.x);
    int xMax = std::max(std::max(v0.x, v1.x), v2.x);
    int yMin = std::min(std::min(v0.y, v1.y), v2.y);
    int yMax = std::max(std::max(v0.y, v1.y), v2.y);

    // Clip bounding box according to nx and ny
    xMin = std::max(std::min(xMin, nx-1), 0);
    xMax = std::min(std::min(xMax, nx-1), nx-1);
    yMin = std::max(std::min(yMin, ny-1), 0);
    yMax = std::min(std::min(yMax, ny-1), ny-1);
    
   //std::cout << xMin << " " << xMax << " " << yMin << " " << yMax << std::endl;

    // Do the overall algorithm from slides
    for(int y=yMin; y<=yMax; y++)
    {
        for(int x=xMin; x<=xMax; x++)
        {
            // Get barycentric coordinates
            // First, compute the area of the triangle formed by the vertices v0, v1, v2
            double area = (v1.y - v2.y) * (v0.x - v2.x) + (v2.x - v1.x) * (v0.y - v2.y);

            double alpha = ((v1.y - v2.y) * (x - v2.x) + (v2.x - v1.x) * (y - v2.y)) / area;
            double beta  = ((v2.y - v0.y) * (x - v2.x) + (v0.x - v2.x) * (y - v2.y)) / area;
            double gamma = 1.0 - alpha - beta;

            //std::cout << alpha << " " << beta << " " << gamma << std::endl;

            // Check if inside triangle
            if(alpha >= 0 && beta >= 0 && gamma >= 0)
            {
                // Should do z-buffering here
                double z = alpha * v0.z + beta * v1.z + gamma * v2.z;  // Check this!!
                //std::cout << "Z value : " << z << std::endl;
                if(depth[x][y] > z)
                {
                    depth[x][y] = z;
                    double r = std::max(0.0, std::min(255.0, alpha * c0->r + beta * c1->r + gamma * c2->r));
                    double g = std::max(0.0, std::min(255.0, alpha * c0->g + beta * c1->g + gamma * c2->g));
                    double b = std::max(0.0, std::min(255.0, alpha * c0->b + beta * c1->b + gamma * c2->b));
                    image[x][y] = Color(r, g, b);
                }
            }
        }
    }
}
// CHECK THIS
bool visible(double den, double num, double &te, double &tl)
{
    // From Lecture Slides:
    if(den > 0.0)
    {
        double t = num / den;
        if(t > tl)
        {
            return false;
        }
        if(t > te)
        {
            te = t;
        }
    }
    else if(den < 0.0)
    {
        double t = num / den;
        if(t < te)
        {
            return false;
        }
        if(t < tl)
        {
            tl = t;
        }
    }
    else
    {
        if(num > 0.0)
        {
            return false;
        }
    }
    return true;
}
// CHECK THIS
bool clippingLiangBarsky(Vec4& v0, Vec4& v1, Color& c0, Color& c1)
{
    bool visibility = false;
    double te = 0.0;
    double tl = 1.0;
    double dx = v1.x - v0.x;
    double dy = v1.y - v0.y;
    double dz = v1.z - v0.z;
    Color dc = Color((c1.r - c0.r), (c1.g - c0.g), (c1.b - c0.b));
    double xMin = -1.0;
    double yMin = -1.0;
    double zMin = -1.0;
    double xMax = 1.0;
    double yMax = 1.0;
    double zMax = 1.0;
    if(visible(dx, xMin - v0.x, te, tl) && visible(-dx, v0.x - xMax, te, tl) && visible(dy, yMin - v0.y, te, tl) && visible(-dy, v0.y - yMax, te, tl) && visible(dz, zMin - v0.z, te, tl) && visible(-dz, v0.z - zMax, te, tl))
    {
        //std::cout << "TL: " << tl << " TE: " << te << std::endl;
        if(tl < 1.0)
        {
            v1.x = v0.x + tl * dx;
            v1.y = v0.y + tl * dy;
            v1.z = v0.z + tl * dz;
            c1.r = c0.r + tl * dc.r;
            c1.g = c0.g + tl * dc.g;
            c1.b = c0.b + tl * dc.b;
        }
        if(te > 0.0)
        {
            v0.x = v0.x + te * dx;
            v0.y = v0.y + te * dy;
            v0.z = v0.z + te * dz;
            c0.r = c0.r + te * dc.r;
            c0.g = c0.g + te * dc.g;
            c0.b = c0.b + te * dc.b;
        }
        visibility = true;
    }
    else
    {
        visibility = false;
    }
    return visibility;
}
// CHECK THIS
void lineRasterization(std::vector<std::vector<Color>> & image, std::vector<std::vector<double> > &depth, Vec4 v0, Vec4 v1, Color c0, Color c1) 
{
    // Case 0 : Slope of line is between 0 and 1
    // Case 1 : Slope of line is between 1 and infinity
    // Case 2 : Slope of line is between -1 and 0
    // Case 3 : Slope of line is between -infinity and -1

    // Calculate slope
    double dx = v1.x - v0.x;
    double dy = v1.y - v0.y;
    if(dx == 0)
    {
        dx = dx + EPSILON;
    }
    double slope = dy / dx;
    //std::cout << "Slope : " << slope << std::endl;
    // Case 0 : Increment x , we are going
    if(slope > 0.0 && slope <= 1.0) // dx > dy both positive, or dy > dx both negative
    {
        // Case 0
        if(v1.x < v0.x)
        {
            std::swap(v0, v1);
            std::swap(c0, c1);
            dx = v1.x - v0.x;
        } 
        int y = v0.y;
        double d = 2 * (v0.y - v1.y) + (v1.x - v0.x);
        Color c = c0;
        Color dc = Color((c1.r - c0.r)/dx, (c1.g - c0.g)/dx, (c1.b - c0.b)/dx);
        for (int x = v0.x; x<=v1.x; x++) {
        // Should do z-buffering here
            double t = (x - v0.x) / dx; // Normalized distance between v0.x and v1.x
            double z = (1.0 - t) * v0.z + t * v1.z;
            if(z < depth[x][y])
            {
                depth[x][y] = z;
                double r = std::max(0.0, std::min(255.0, c.r));
                double g = std::max(0.0, std::min(255.0, c.g));
                double b = std::max(0.0, std::min(255.0, c.b));
                image[x][y] = Color(r, g, b);
            }
            if(d < 0.0)
            {
                y += 1;
                d += (v0.y - v1.y) + (v1.x - v0.x);
            }
            else
            {
                d += (v0.y - v1.y);
            }
            c.r += dc.r;
            c.g += dc.g;
            c.b += dc.b;
        }
    }
    else if(slope > 1.0)
    {
        // Case 1
        if(v1.y < v0.y)
        {
            std::swap(v0, v1);
            std::swap(c0, c1);
            dy = v1.y - v0.y;
        }
        int x = v0.x;
        double d = 2.0 * (v0.x - v1.x) + (v1.y - v0.y);
        Color c = c0;
        Color dc = Color((c1.r - c0.r)/dy, (c1.g - c0.g)/dy, (c1.b - c0.b)/dy);
        for (int y = v0.y; y<=v1.y; y++) {
        // Should do z-buffering here
            double t = (y - v0.y) / dy; // Normalized distance between v0.y and v1.y
            double z = (1.0 - t) * v0.z + t * v1.z;
            if(z < depth[x][y])
            {
                depth[x][y] = z;
                double r = std::max(0.0, std::min(255.0, c.r));
                double g = std::max(0.0, std::min(255.0, c.g));
                double b = std::max(0.0, std::min(255.0, c.b));
                image[x][y] = Color(r, g, b);
            }
            if(d < 0.0)
            {
                x += 1;
                d += (v0.x - v1.x) + (v1.y - v0.y);
            }
            else
            {
                d += (v0.x - v1.x);
            }
            c.r += dc.r;
            c.g += dc.g;
            c.b += dc.b;
        }
    }
    else if(slope >= -1.0 && slope <= 0.0)
    {
        // Case 2
        if(v1.x < v0.x)
        {
            std::swap(v0, v1);
            std::swap(c0, c1);
            dx = v1.x - v0.x;
        }
        int y = v0.y;
        double d = 2.0 * (v0.y - v1.y) - (v1.x - v0.x);
        Color c = c0;
        Color dc = Color((c1.r - c0.r)/dx, (c1.g - c0.g)/dx, (c1.b - c0.b)/dx);
        for (int x = v0.x; x<=v1.x; x++) {
        // Should do z-buffering here
            double t = (x - v0.x) / dx; // Normalized distance between v0.x and v1.x
            double z = (1.0 - t) * v0.z + t * v1.z;
            if(z < depth[x][y])
            {
                depth[x][y] = z;
                double r = std::max(0.0, std::min(255.0, c.r));
                double g = std::max(0.0, std::min(255.0, c.g));
                double b = std::max(0.0, std::min(255.0, c.b));
                image[x][y] = Color(r, g, b);
            }
            if(d > 0.0)
            {
                y -= 1;
                d += (v0.y - v1.y) - (v1.x - v0.x);
            }
            else
            {
                d += (v0.y - v1.y);
            }
            c.r += dc.r;
            c.g += dc.g;
            c.b += dc.b;
        }
    }
    else if(slope < -1.0)
    {
        // Case 3
        if(v1.y < v0.y)
        {
            std::swap(v0, v1);
            std::swap(c0, c1);
            dy = v1.y - v0.y;
        }
        int x = v0.x;
        double d = 2.0 * (v0.x - v1.x) - (v1.y - v0.y);
        Color c = c0;
        Color dc = Color((c1.r - c0.r)/dy, (c1.g - c0.g)/dy, (c1.b - c0.b)/dy);
        for (int y = v0.y; y<=v1.y; y++) {
        // Should do z-buffering here
            double t = (y - v0.y) / dy; // Normalized distance between v0.y and v1.y
            double z = (1.0 - t) * v0.z + t * v1.z;
            if(z < depth[x][y])
            {
                depth[x][y] = z;
                double r = std::max(0.0, std::min(255.0, c.r));
                double g = std::max(0.0, std::min(255.0, c.g));
                double b = std::max(0.0, std::min(255.0, c.b));
                image[x][y] = Color(r, g, b);
            }
            if(d > 0.0)
            {
                x -= 1;
                d += (v0.x - v1.x) - (v1.y - v0.y);
            }
            else
            {
                d += (v0.x - v1.x);
            }
            c.r += dc.r;
            c.g += dc.g;
            c.b += dc.b;
        }
    }
}

