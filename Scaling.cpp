#include "Scaling.h"
#include <iostream>
#include <iomanip>
#include "Matrix4.h"

using namespace std;

Scaling::Scaling() {}

Scaling::Scaling(int scalingId, double sx, double sy, double sz)
{
    this->scalingId = scalingId;
    this->sx = sx;
    this->sy = sy;
    this->sz = sz;
}

Matrix4 Scaling::getScalingMatrix() {
    double scalingMatrix[4][4] = {
        { this->sx, 0, 0, 0 },
        { 0, this->sy, 0, 0 },
        { 0, 0, this->sz, 0},
        { 0, 0, 0, 1 }
    };

    return Matrix4(scalingMatrix);
}

ostream &operator<<(ostream &os, const Scaling &s)
{
    os << fixed << setprecision(3) << "Scaling " << s.scalingId << " => [" << s.sx << ", " << s.sy << ", " << s.sz << "]";

    return os;
}
