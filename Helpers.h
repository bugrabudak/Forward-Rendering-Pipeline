#ifndef __HELPERS_H__
#define __HELPERS_H__

#define ABS(a) ((a) > 0 ? (a) : -1 * (a))
#define EPSILON 0.000000001

#include <vector>
#include "Matrix4.h"
#include "Vec3.h"
#include "Vec4.h"
#include "Camera.h"
#include "Color.h"
#include "Mesh.h"
#include "Translation.h"
#include "Scaling.h"
#include "Rotation.h"

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

Matrix4 getViewportTransformationMatrix(Camera*);

Matrix4 getCameraTransformationMatrix(Camera*);

Matrix4 getOrthographicProjectionTransformationMatrix(Camera * camera);

Matrix4 getPerspectiveProjectionTransformationMatrix(Camera * camera);

bool visible(double den, double num, double &tE, double &tL);

bool clip(double xMax, double xMin, double yMax, double yMin, double zMax, double zMin, Vec4 &vec1, Vec4 &vec2, Color* color1, Color* color2);

Matrix4 getModelingTransformationMatrix(Camera* camera, Mesh* mesh, vector<Translation*>& translations, vector<Scaling*>& scalings, vector<Rotation*>& rotations);

Matrix4 getCumulativeTransformations(Matrix4&, Matrix4&, Matrix4&);

bool isCulled(Vec4* transformedArray);

void rasterize(vector< vector<Color> >& image, Vec4 first, Vec4 second, Color firstColor, Color secondColor);

void lineRasterization(vector< vector<Color> >& image, Vec4 first, Vec4 second, Color firstColor, Color secondColor);

#endif