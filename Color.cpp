#include "Color.h"
#include <iostream>
#include <iomanip>

using namespace std;

Color::Color() {}

Color::Color(double r, double g, double b)
{
    this->r = r;
    this->g = g;
    this->b = b;
}

Color::Color(const Color &other)
{
    this->r = other.r;
    this->g = other.g;
    this->b = other.b;
}

Color Color::operator-(Color& second) {
    Color result;
    result.r = this->r - second.r;
    result.g = this->g - second.g;
    result.b = this->b - second.b;
    return result;
}

Color Color::operator/(double second) {
    Color result;
    result.r = this->r / second;
    result.g = this->g / second;
    result.b = this->b / second;
    return result;
}

Color Color::operator+(Color& second) {
    Color result;
    result.r = this->r + second.r;
    result.g = this->g + second.g;
    result.b = this->b + second.b;
    return result;
}


ostream& operator<<(ostream& os, const Color& c)
{
    os << fixed << setprecision(0) << "rgb(" << c.r << ", " << c.g << ", " << c.b << ")";
    return os;
}
