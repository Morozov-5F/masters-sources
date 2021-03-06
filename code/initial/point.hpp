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

    double GetX() const;
    double GetY() const;
    double GetZ() const;

    void SetX(double val);
    void SetY(double val);
    void SetZ(double val);

    void Normalize();

    double GetMagnitude() const;
    double GetSquaredMagnitude() const;

    double GetPolar() const;
    double GetIncl() const;

    double Dot(const Point& point) const;

    double GetAngle(const Point& point) const;

    Point GetRotated(const Point& euler_angles) const;

    Point operator+(const Point& point) const;
    Point operator-(const Point& point) const;
    Point operator*(double scalar) const;

    Point& operator+=(const Point& point);
    Point& operator-=(const Point& point);
    Point& operator*=(double scalar);

    friend std::istream& operator>>(std::istream& is, Point& point);
    friend std::ostream& operator<<(std::ostream& os, const Point& point);
};

#endif // POINT_H
