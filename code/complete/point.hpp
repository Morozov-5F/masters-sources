#ifndef POINT_H
#define POINT_H

#include <fstream>

class Point
{
private:
    double x;
    double y;
    double z;

public:
    Point() = default;
    Point(double x, double y, double z);

    static double GetDistance(const Point& point1, const Point& point2);
    static double GetSquaredDistance(const Point& point1,
                                     const Point& point2);
    static double GetCosBetween(const Point& p1, const Point& p2);

    double GetX() const;
    double GetY() const;
    double GetZ() const;

    void SetX(double val);
    void SetY(double val);
    void SetZ(double val);

    void Normalize();
    Point Normalized() const;

    double GetMagnitude() const;
    double GetSquaredMagnitude() const;

    double GetPolar() const;    // TODO: Specify coordinate system type (left/right)
    double GetAzimuth() const;  // TODO: Specify coordinate system type (left/right)

    double Dot(const Point& point) const;
    Point Cross(const Point& point) const;

    double GetAngle(const Point& point) const;

    Point GetRotated(const Point& euler_angles) const; // TODO: Specify the coordinate system type

    Point operator-() const;

    Point operator+(const Point& point) const;
    Point operator-(const Point& point) const;
    Point operator*(double scalar) const;
    Point operator/(double scalar) const;

    Point& operator+=(const Point& point);
    Point& operator-=(const Point& point);
    Point& operator*=(double scalar);

    friend std::istream& operator>>(std::istream& is, Point& point);
    friend std::ostream& operator<<(std::ostream& os, const Point& point);

    friend Point operator*(double scalar, const Point& point);
};

namespace geom_constants {
    static const Point OX = Point(1, 0, 0);
    static const Point OY = Point(0, 1, 0);
    static const Point OZ = Point(0, 0, 1);
}

#endif // POINT_H
