#include "Rotation.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include "Vec3.h"

using namespace std;

Vec3 crossProduct(Vec3 a, Vec3 b)
{
    Vec3 result;

    result.x = a.y * b.z - b.y * a.z;
    result.y = b.x * a.z - a.x * b.z;
    result.z = a.x * b.y - b.x * a.y;

    return result;
}

Vec3 normalize(Vec3 v)
{
    Vec3 result;
    double d;

    d = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);

    result.x = v.x / d;
    result.y = v.y / d;
    result.z = v.z / d;

    return result;
}


Rotation::Rotation() {}

Rotation::Rotation(int rotationId, double angle, double x, double y, double z)
{
    this->rotationId = rotationId;
    this->angle = angle;
    this->ux = x;
    this->uy = y;
    this->uz = z;
}

Matrix4* Rotation::getRotationMatrices() {
    Matrix4* arr = new Matrix4[3];
    double smallest = min(abs(this->ux), min(abs(this->uy), abs(this->uz)));
    Vec3 v;
    Vec3 u = Vec3(this->ux, this->uy, this->uz, 0);

    if (smallest == abs(this->ux)) {
        v = Vec3(0, -(this->uz), this->uy, 0);
    } else if (smallest == abs(this->uy)) {
        v =  Vec3(-(this->uz), 0, this->ux, 0);
    } else {
        v = Vec3(-(this->uy), this->ux, 0, 0);
    }

    Vec3 w = crossProduct(u, v);
    v = normalize(v);
    w = normalize(w);

    double transformMatrix[4][4] = {
        { u.x, u.y, u.z, 0 },
        { v.x, v.y, v.z, 0 },
        { w.x, w.y, w.z, 0 },
        { 0, 0, 0, 1 }
    };

    double transformMatrixT[4][4] = {
        { u.x, v.x, w.x, 0 },
        { u.y, v.y, w.y, 0 },
        { u.z, v.z, w.z, 0 },
        { 0, 0, 0, 1 }
    };

    double rotationMatrix[4][4] = {
        { 1, 0, 0, 0 },
        { 0, cos(this->angle * M_PI / 180), -sin(this->angle * M_PI / 180), 0 },
        { 0, sin(this->angle * M_PI / 180), cos(this->angle * M_PI / 180), 0 },
        { 0, 0, 0, 1 }
    };

    arr[0] = Matrix4(transformMatrix);
    arr[1] = Matrix4(rotationMatrix);
    arr[2] = Matrix4(transformMatrixT);

    return arr;  
}

ostream &operator<<(ostream &os, const Rotation &r)
{
    os << fixed << setprecision(3) << "Rotation " << r.rotationId << " => [angle=" << r.angle << ", " << r.ux << ", " << r.uy << ", " << r.uz << "]";

    return os;
}