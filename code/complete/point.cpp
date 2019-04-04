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
    auto normalized = Normalized();

    x /= normalized.x;
    y /= normalized.y;
    z /= normalized.z;
}

Point Point::Normalized() const
{
    double mag = GetMagnitude();

    if (!mag) return {0, 0, 0};

    return {x / mag, y / mag, z / mag};
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

double Point::GetCosBetween(const Point& p1, const Point& p2)
{
    return p1.GetAngle(p2);
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

double Point::GetAzimuth() const
{
    return std::asin(z / GetMagnitude());
}

double Point::Dot(const Point& point) const
{
    return x * point.x + y * point.y + z * point.z;
}

Point Point::Cross(const Point& point) const
{
    double nx = y * point.z - z * point.y;
    double ny = z * point.x - x * point.z;
    double nz = x * point.y - y * point.x;

    return {nx, ny, nz};
}

double Point::GetAngle(const Point& point) const
{
    double phi_1 = GetPolar(), phi_2 = point.GetPolar();
    double theta_1 = GetAzimuth(), theta_2 = point.GetAzimuth();

    Point a(std::cos(phi_1) * std::cos(theta_1),
            std::sin(phi_1) * std::cos(theta_1),
            std::sin(theta_1));

    Point b(std::cos(phi_2) * std::cos(theta_2),
            std::sin(phi_2) * std::cos(theta_2),
            std::sin(theta_2));

    return a.Dot(b);
}

Point Point::GetRotated(const Point& euler_angles) const
{
    double cx = std::cos(euler_angles.GetX()), sx = std::sin(euler_angles.GetX());
    double cy = std::cos(euler_angles.GetY()), sy = std::sin(euler_angles.GetY());
    double cz = std::cos(euler_angles.GetZ()), sz = std::sin(euler_angles.GetZ());

    return Point(x * cz * cy - y * cy * sz + z * sy,
                 x * (cx * sz + cz * sx * sy) + y * (cz * cx - sz * sx * sy) - z * cy * sx,
                 x * (sz * sx - cz * cx * sy) + y * (cz * sx + cx * sz * sy) + z * cy * cx);
}

Point Point::operator-() const {
    return {-x, -y, -z};
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

Point Point::operator/(double scalar) const
{
    if (scalar == 0) return {NAN, NAN, NAN};

    return { x / scalar, y / scalar, z / scalar };
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

std::istream& operator>>(std::istream& is, Point& point)
{
    is >> point.x >> point.y >> point.z;
    return is;
}

std::ostream& operator<<(std::ostream& os, const Point& point)
{
    os << "(" << point.x << " " << point.y << " " << point.z << ")";
    return os;
}

Point operator+(double scalar, const Point& point)
{
    return point * scalar;
}
