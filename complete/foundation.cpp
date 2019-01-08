#include "foundation.hpp"

#include <iostream>

#include <cmath>

namespace GC = geom_constants;

Foundation::Foundation(double R, double angle1, double angle2, double p1, double p2) : r(R)
{
    auto angle3 = 180. - angle1 - angle2;

    if (angle3 <= 0.)
    {
        std::cerr << "Icorrect angles: should be less than 180 degrees" << std::endl;
        return;
    }

    if (angle3 >= 90.)
    {
        this->alpha = std::max(angle1, angle2);
        this->beta = angle3;
        this->gamma = std::min(angle1, angle2);
    }
    else
    {
        this->alpha = std::max(angle3, std::min(angle1, angle2));
        this->beta = std::max(angle1, angle2);
        this->gamma = std::min(angle3, std::min(angle1, angle2));
    }

    this->alpha *= M_PI / 180.;
    this->beta *= M_PI / 180.;
    this->gamma *= M_PI / 180.;

    // Not optimized, needed for debugging
    double a = 2 * R * std::sin(this->alpha);
    double b = 2 * R * std::sin(this->beta);
    double c = 2 * R * std::sin(this->gamma);

    //
    // Geometry of our triangle
    //
    //         M3
    //        *.
    //       .    .  a
    //    b .         .
    //     .             .
    //    *. . . . . . . . *
    //   M1        c        M2
    //
    // Alpha is M3M1M2
    // Beta is M1M2M3
    // Gamma is M1M3M2
    //

    // Find coordinates keeping in mind that z coordinates are set using
    // p1 and p2 coefficients in a following way:
    // - M2(x2, y2, p1 * R)
    // - M3(x3, y3, p2 * R)

    auto k1 = p1 * R, k2 = p2 * R;

    this->m1 = Point();
    this->m2 = Point(std::sqrt(c * c - k1 * k1), 0, k1);

    auto x = (c * c + b * b - a * a - 2 * k1 * k2) / (2 * std::sqrt(c * c - k1 * k1));
    this->m3 = Point(x, std::sqrt(b * b - k2 * k2 - x * x), k2);
    // Kind of a self-check
    this->d12 = Point::GetDistance(this->m1, this->m2);
    this->d23 = Point::GetDistance(this->m2, this->m3);
    this->d13 = Point::GetDistance(this->m1, this->m3);

    // Find transform components for local coordinate system
    this->center = (this->m1 + this->m2 + this->m3) / 3;
    auto n = -(this->m2 - this->m1).Cross(this->m3 - this->m1);
    if (n.Dot(GC::OZ) < 0)
    {
        n = -n;
    }

    this->e1 = (this->m1 - this->center).Normalized();
    this->e3 = n.Normalized();
    this->e2 = -this->e1.Cross(this->e3).Normalized();

    this->m1_local = this->GlobalToLocal(this->m1);
    this->m2_local = this->GlobalToLocal(this->m2);
    this->m3_local = this->GlobalToLocal(this->m3);

    this->PrintInfo();

#ifndef NDEBUG
    std::cout << "Edges (initial):" << std::endl;
    std::cout << "  d12 (M1.M2, c): " << c << "m" << std::endl;
    std::cout << "  d13 (M1.M3, b): " << b << "m" << std::endl;
    std::cout << "  d23 (M2.M3, a): " << a << "m" << std::endl;
#endif
}

Point Foundation::LocalToGlobal(const Point& point) const
{
    // TODO: Optimize calls
    double x = point.GetX() * e1.Dot(GC::OX) + point.GetY() * e2.Dot(GC::OX) + point.GetZ() * e3.Dot(GC::OX);
    double y = point.GetX() * e1.Dot(GC::OY) + point.GetY() * e2.Dot(GC::OY) + point.GetZ() * e3.Dot(GC::OY);
    double z = point.GetX() * e1.Dot(GC::OZ) + point.GetY() * e2.Dot(GC::OZ) + point.GetZ() * e3.Dot(GC::OZ);

    return center + Point(x, y, z);
}

Point Foundation::GlobalToLocal(const Point& point) const
{
    auto temp = point - center;

    double x = temp.GetX() * GC::OX.Dot(e1) + temp.GetY() * GC::OY.Dot(e1) + temp.GetZ() * GC::OZ.Dot(e1);
    double y = temp.GetX() * GC::OX.Dot(e2) + temp.GetY() * GC::OY.Dot(e2) + temp.GetZ() * GC::OZ.Dot(e2);
    double z = temp.GetX() * GC::OX.Dot(e3) + temp.GetY() * GC::OY.Dot(e3) + temp.GetZ() * GC::OZ.Dot(e3);

    return {x, y, z};
}

void Foundation::PrintInfo() const
{
    std::cout << "Foundation information:" << std::endl;
    std::cout << "Local center:" << std::endl;
    std::cout << "  O\": " << center << std::endl;
    std::cout << "Coordinates (global):" << std::endl;
    std::cout << "  M1: " << m1 << std::endl;
    std::cout << "  M2: " << m2 << std::endl;
    std::cout << "  M3: " << m3 << std::endl;
#ifndef NDEBUG
    std::cout << "Coordinates (global, deduced from local):" << std::endl;
    std::cout << "  M1: " << LocalToGlobal(m1_local) << std::endl;
    std::cout << "  M2: " << LocalToGlobal(m2_local) << std::endl;
    std::cout << "  M3: " << LocalToGlobal(m3_local) << std::endl;
    std::cout << "Coordinates (local, deduced from global):" << std::endl;
    std::cout << "  M1: " << GlobalToLocal(m1) << std::endl;
    std::cout << "  M2: " << GlobalToLocal(m2) << std::endl;
    std::cout << "  M3: " << GlobalToLocal(m3) << std::endl;
#endif
    std::cout << "Coordinates (local):" << std::endl;
    std::cout << "  M1: " << m1_local << std::endl;
    std::cout << "  M2: " << m2_local << std::endl;
    std::cout << "  M3: " << m3_local << std::endl;
    std::cout << "Radius:" << std::endl;
    std::cout << "  R: " << r << " m" << std::endl;
    std::cout << "Angles:" << std::endl;
    std::cout << "  M2.M1.M3: " << alpha * 180 / M_PI << " deg" << std::endl;
    std::cout << "  M1.M2.M3: " << beta * 180 / M_PI << " deg" << std::endl;
    std::cout << "  M1.M3.M2: " << gamma * 180 / M_PI << " deg" << std::endl;
    std::cout << "Edges:" << std::endl;
    std::cout << "  d12 (M1.M2): " << d12 << " m" << std::endl;
    std::cout << "  d13 (M1.M3): " << d13 << " m" << std::endl;
    std::cout << "  d23 (M2.M3): " << d23 << " m" << std::endl;
#ifndef NDEBUG
    std::cout << "Edges (from local):" << std::endl;
    std::cout << "  d12 (M1.M2): " << Point::GetDistance(m1_local, m2_local) << " m" << std::endl;
    std::cout << "  d13 (M1.M3): " << Point::GetDistance(m1_local, m3_local) << " m" << std::endl;
    std::cout << "  d23 (M2.M3): " << Point::GetDistance(m2_local, m3_local)  << " m"<< std::endl;
#endif
}
