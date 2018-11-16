#include "point.hpp"
#include <cmath>

Point::Point(double x, double y, double z): x(x), y(y), z(z)
{}

double Point::GetX() const
{
    return x;
}

double Point::GetY() const
{
    return y;
}

double Point::GetZ() const
{
    return z;
}

void Point::Normalize()
{
    double mag = GetMagnitude();

    if (!mag) return;

    x /= mag;
    y /= mag;
    z /= mag;
}

void Point::SetX(double val)
{
    x = val;
}

void Point::SetY(double val)
{
    y = val;
}

void Point::SetZ(double val)
{
    z = val;
}

double Point::GetDistance(const Point& point1, const Point& point2)
{
    return (point1 - point2).GetMagnitude();
}

double Point::GetSquaredDistance(const Point& point1,
                                        const Point& point2)
{
    return (point1 - point2).GetSquaredMagnitude();
}

double Point::GetSquaredMagnitude() const
{
    return x * x + y * y + z * z;
}

double Point::GetMagnitude() const
{
    return std::sqrt(GetSquaredMagnitude());
}

double Point::GetPolar() const
{
    return std::atan2(y, x);
}

double Point::GetIncl() const
{
    return std::asin(z / GetMagnitude());
}

double Point::Dot(const Point& point) const
{
    return x * point.x + y * point.y + z * point.z;
}

double Point::GetAngle(const Point& point) const
{
    double phi_1 = GetPolar(), phi_2 = point.GetPolar();
    double theta_1 = GetIncl(), theta_2 = point.GetIncl();

    Point a(std::cos(phi_1) * std::cos(theta_1),
            std::sin(phi_1) * std::cos(theta_1),
            std::sin(theta_1));

    Point b(std::cos(phi_2) * std::cos(theta_2),
            std::sin(phi_2) * std::cos(theta_2),
            std::sin(theta_2));

    return a.Dot(b);
}

Point Point::operator+(const Point& point) const
{
    return { x + point.x, y + point.y, z + point.z };
}

Point Point::operator-(const Point& point) const
{
    return { x - point.x, y - point.y, z - point.z };
}

Point Point::operator*(double scalar) const
{
    return { x * scalar, y * scalar, z * scalar };
}

Point& Point::operator+=(const Point& point)
{
    this->x += point.x;
    this->y += point.y;
    this->y += point.z;

    return *this;
}

Point& Point::operator-=(const Point& point)
{
    this->x -= point.x;
    this->y -= point.y;
    this->y -= point.z;

    return *this;
}

Point& Point::operator*=(double scalar)
{
    this->x *= scalar;
    this->y *= scalar;
    this->y *= scalar;

    return *this;
}
