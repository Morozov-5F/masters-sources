#pragma once

#include "point.hpp"

class Foundation
{
private:
    Point m1, m2, m3;                     /// Vertices
    Point m1_local, m2_local, m3_local;   /// Vertices in the global coordinate system

    Point e1, e2, e3;           // Local coordinate system basis
    Point center;               // Coordinate system center

    double d12, d13, d23;       /// Side lengths (for optimization?)

    double r;                   /// Radius of the circumscribed circle
    double alpha, beta, gamma;  /// Triangle angles

public:
    Foundation(double R, double angle1, double angle2, double m1, double m2);

    Point GlobalToLocal(const Point& point) const;
    Point LocalToGlobal(const Point& point) const;

    void PrintInfo() const;
};
