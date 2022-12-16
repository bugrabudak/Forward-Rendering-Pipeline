#include <iostream>
#include <cmath>
#include "Helpers.h"
#include "Matrix4.h"
#include "Vec3.h"
#include "Vec4.h"
#include "Camera.h"
#include "Color.h"

using namespace std;

/*
 * Calculate cross product of vec3 a, vec3 b and return resulting vec3.
 */
Vec3 crossProductVec3(Vec3 a, Vec3 b)
{
    Vec3 result;

    result.x = a.y * b.z - b.y * a.z;
    result.y = b.x * a.z - a.x * b.z;
    result.z = a.x * b.y - b.x * a.y;

    return result;
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
    Vec3 result;
    double d;

    d = magnitudeOfVec3(v);
    result.x = v.x / d;
    result.y = v.y / d;
    result.z = v.z / d;

    return result;
}

/*
 * Return -v (inverse of vec3 v)
 */
Vec3 inverseVec3(Vec3 v)
{
    Vec3 result;
    result.x = -v.x;
    result.y = -v.y;
    result.z = -v.z;

    return result;
}

/*
 * Add vec3 a to vec3 b and return resulting vec3 (a+b).
 */
Vec3 addVec3(Vec3 a, Vec3 b)
{
    Vec3 result;
    result.x = a.x + b.x;
    result.y = a.y + b.y;
    result.z = a.z + b.z;

    return result;
}

/*
 * Subtract vec3 b from vec3 a and return resulting vec3 (a-b).
 */
Vec3 subtractVec3(Vec3 a, Vec3 b)
{
    Vec3 result;
    result.x = a.x - b.x;
    result.y = a.y - b.y;
    result.z = a.z - b.z;

    return result;
}

/*
 * Multiply each element of vec3 with scalar.
 */
Vec3 multiplyVec3WithScalar(Vec3 v, double c)
{
    Vec3 result;
    result.x = v.x * c;
    result.y = v.y * c;
    result.z = v.z * c;

    return result;
}

/*
 * Prints elements in a vec3. Can be used for debugging purposes.
 */
void printVec3(Vec3 v)
{
    cout << "(" << v.x << "," << v.y << "," << v.z << ")" << endl;
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
                result.val[i][j] = 1.0;
            }
            else
            {
                result.val[i][j] = 0.0;
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
                total += m1.val[i][k] * m2.val[k][j];
            }

            result.val[i][j] = total;
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
            total += m.val[i][j] * v.getElementAt(j);
        }
        values[i] = total;
    }

    return Vec4(values[0], values[1], values[2], values[3], v.colorId);
}

Matrix4 getViewportTransformationMatrix(Camera * camera) {
    double viewportArray[4][4] = {
        { camera->horRes / 2.0, 0, 0, ((camera->horRes - 1) / 2.0) + camera->left },
        { 0, camera->verRes / 2.0, 0, ((camera->verRes - 1) / 2.0) - camera->bottom },
        { 0, 0, 0.5, 0.5 },
        { 0, 0, 0, 1 }
    };

    Matrix4 viewportTransformationMatrix = Matrix4(viewportArray);
    return viewportTransformationMatrix;
}

Matrix4 getCameraTransformationMatrix(Camera * camera) {
    double firstArray[4][4] = {{1, 0, 0, -(camera->pos.x)},
                    {0, 1, 0, -(camera->pos.y)},
                    {0, 0, 1, -(camera->pos.z)},
                    {0, 0, 0, 1}};
    double secondArray[4][4] = {{camera->u.x, camera->u.y, camera->u.z, 0},
                        {camera->v.x, camera->v.y, camera->v.z, 0},
                        {camera->w.x, camera->w.y, camera->w.z, 0},
                        {0, 0, 0, 1}};
    Matrix4 first = Matrix4(firstArray);
    Matrix4 second = Matrix4(secondArray);
    return multiplyMatrixWithMatrix(first, second);
}

Matrix4 getPerspectiveProjectionTransformationMatrix(Camera * camera) {
    double m1 = (2*camera->near) / (camera->right - camera->left);
    double m3 = (camera->right + camera->left) / (camera->right - camera->left);
    double m6 = (2*camera->near) / (camera->top - camera->bottom); 
    double m7 = (camera->top + camera->bottom) / (camera->top - camera->bottom);
    double m10 =  -((camera->far + camera->near) / (camera->far - camera->near));
    double m11 = -((2*camera->far*camera->near) / (camera->far - camera->near));

    double perspectiveArray[4][4] = {{m1, 0, m3, 0},
                                    {0, m6, m7, 0},
                                    {0, 0, m10, m11},
                                    {0, 0, -1, 0}};
    Matrix4 perspectiveMatrix = Matrix4(perspectiveArray);
    return perspectiveMatrix;
}

Matrix4 getOrthographicProjectionTransformationMatrix(Camera * camera) {
    double m1 = 2 / (camera->right - camera->left);
    double m4 = -((camera->right + camera->left) / (camera->right - camera->left));
    double m6 = 2/(camera->top - camera->bottom);
    double m8 = -((camera->top + camera->bottom) / (camera->top - camera->bottom));
    double m11 = -(2/(camera->far - camera->near));
    double m12 = -((camera->far + camera->near) / (camera->far - camera->near));

    double orthographicArray[4][4] = {{m1, 0, 0, m4},
                                        {0, m6, 0, m8},
                                        {0, 0, m11,m12},
                                        {0, 0, 0, 1}};

    Matrix4 orthographicMatrix = Matrix4(orthographicArray);
    return orthographicMatrix;
}

// visible code from the slides
bool visible(double den, double num, double &tE, double &tL){
    double t = num/den;

    if (den > 0) {        
        if (t > tL){
            return false;
        }
        if (t > tE) {
            tE = t;
        }
    } else if (den < 0) {    
        if (t < tE) {
            return false;
        }
        if (t < tL) {
            tL = t;
        }
    }
    else if (num > 0) {
        return false;
    }

    return true;
}

bool clip(double xMax, double xMin, double yMax, double yMin, double zMax, double zMin, 
          Vec4 &vec1, Vec4 &vec2, Color &color1, Color &color2){
    double dX, dY, dZ;
    double tE = 0;
    double tL = 1;
    bool isVisible = false;

    dX = vec2.x - vec1.x;
    dY = vec2.y - vec1.y;
    dZ = vec2.z - vec1.z;

    Color dColor = Color(color2.r - color1.r,color2.g - color1.g, color2.b - color1.b);

    if (visible(dX, xMin - vec1.x, tE, tL) && visible(-dX, vec1.x - xMax, tE, tL) &&
        visible(dY, yMin - vec1.y, tE, tL) && visible(-dY, vec1.y - yMax, tE, tL) && 
        visible(dZ, zMin - vec1.z, tE, tL) && visible(-dZ, vec1.z - zMax, tE, tL)){

        if (tL < 1){
            vec2.x = vec1.x + dX * tL;
            vec2.y = vec1.y + dY * tL;
            vec2.z = vec1.z + dZ * tL;
            double newR = color1.r + (dColor.r * tL);
            double newG = color1.g + (dColor.g * tL);
            double newB = color1.b + (dColor.b * tL);
            color2.r = newR;
            color2.g = newG;
            color2.b= newB;
        }
        if (tE > 0){
            vec1.x = vec1.x + dX * tE;
            vec1.y = vec1.y + dY * tE;
            vec1.z = vec1.z + dZ * tE;

            double newR = color1.r + (dColor.r * tE);
            double newG = color1.g + (dColor.g * tE);
            double newB = color1.b + (dColor.b * tE);
            color1.r = newR;
            color1.g = newG;
            color1.b= newB;
        }
        isVisible = true;
    }

    return isVisible;
}