#include <iostream>
#include <cmath>
#include <vector>
#include <random>

#include <algorithm>

#include "point.hpp"

double get_radius(double d12, double d13, double d23)
{
    double p = (d12 + d13 + d23) / 2.0;

    return (d12 * d13 * d23) / (4 * std::sqrt(p * (p - d12) * (p - d13) * (p - d23)));
}

void read_data(Point M[4], double &lambda, double &sigma)
{
    // TODO: read data from the file
    double phi = 35, psi = 270, theta = 5;

    M[0] = Point( 0.5, -0.5, 2.0);
    M[1] = Point( 2.0,  1.0, 0.0);
    M[2] = Point(-2.0,  1.0, 0.0);
    M[3] = Point( 1.0, -2.0, 0.0);


    lambda = 0.0;
    sigma = 0.1;
}

Point find_coordinates(Point M[], double l[])
{
    double x1 = M[1].GetX(), y1 = M[1].GetY();
    double x2 = M[2].GetX(), y2 = M[2].GetY();
    double x3 = M[3].GetX(), y3 = M[3].GetY();

    double P = (x1 - x2) * (y1 - y3) - (x1 - x3) * (y1 - y2);
    double R = (x1 - x2) * (l[3] - l[1] + (x1 * x1 - x3 * x3) + (y1 * y1 - y3 * y3)) -
               (x1 - x3) * (l[2] - l[1] + (x1 * x1 - x2 * x2) + (y1 * y1 - y2 * y2));
    double Q = (y1 - y3) * (l[2] - l[1] + (x1 * x1 - x2 * x2) + (y1 * y1 - y2 * y2)) -
               (y1 - y2) * (l[3] - l[1] + (x1 * x1 - x3 * x3) + (y1 * y1 - y3 * y3));

    double x = Q / 2 / P;
    double y = R / 2 / P;
    double z = std::sqrt(l[1] - (x - x1) * (x - x1) - (y - y1) * (y - y1));

    return { x, y, z };
}

double solve(Point M[5], double &err_l, double &err_c)
{
    // foundation sides
    double d12 = Point::GetSquaredDistance(M[1], M[2]);
    double d23 = Point::GetSquaredDistance(M[2], M[3]);
    double d31 = Point::GetSquaredDistance(M[3], M[1]);

#if 0
    Point add = M[1] + (M[2] - M[1]) * (0.5 + lambda * radius);
    add += Point(-add.GetY(), add.GetX(), add.GetZ()) * sigma * radius;

    M[4] = add;
#endif

    double d14 = Point::GetSquaredDistance(M[1], M[4]);
    double d42 = Point::GetSquaredDistance(M[4], M[2]);

    Point dist[] = { Point(),
                     M[1] - M[0],
                     M[2] - M[0],
                     M[3] - M[0],
                     M[4] - M[0]};

    // sides to top
    double l1 = dist[1].GetSquaredMagnitude();
    double l2 = dist[2].GetSquaredMagnitude();
    double l3 = dist[3].GetSquaredMagnitude();
    double l4 = dist[4].GetSquaredMagnitude();

    // cosine of angles
    double cA12_ctrl = (l1 + l2 - d12) / (2 * std::sqrt(l1) * std::sqrt(l2));
    double cA23_ctrl = (l2 + l3 - d23) / (2 * std::sqrt(l2) * std::sqrt(l3));
    double cA31_ctrl = (l3 + l1 - d31) / (2 * std::sqrt(l3) * std::sqrt(l1));
    double cA14_ctrl = (l1 + l4 - d14) / (2 * std::sqrt(l1) * std::sqrt(l4));
    double cA42_ctrl = (l4 + l2 - d42) / (2 * std::sqrt(l4) * std::sqrt(l2));

    double cA12 = dist[1].GetAngle(dist[2]);
    double cA23 = dist[2].GetAngle(dist[3]);
    double cA31 = dist[3].GetAngle(dist[1]);
    double cA14 = dist[1].GetAngle(dist[4]);
    double cA42 = dist[4].GetAngle(dist[2]);

    // sine of 4-point angles
    double sA14 = std::sqrt(1 - cA14 * cA14);
    double sA42 = std::sqrt(1 - cA42 * cA42);

    // we're trying find sides to the apex
    double l1_, l2_, l3_, k = std::sqrt(d42 / d14) * sA14 / sA42;
    l1_ = d12 / (1 + k * k - 2 * k * cA12);
    l2_ = k * k * l1;

    double D1 = l1_ * cA31 * cA31 - l1_ + d31;
    double D2 = l2_ * cA23 * cA23 - l2_ + d23;

    double x11 = std::sqrt(l1_) * cA31 + std::sqrt(D1);
    double x12 = std::sqrt(l1_) * cA31 - std::sqrt(D1);

    double x21 = std::sqrt(l2_) * cA23 + std::sqrt(D2);
    double x22 = std::sqrt(l2_) * cA23 - std::sqrt(D2);

    double sol1 = std::max(x11, x12);
    double sol2 = std::max(x21, x22);
    l3_ = (sol1 + sol2) / 2;
    l3_ *= l3_;

    double l_new[] = { 0, l1_, l2_, l3_ };
    auto coord = find_coordinates(M, l_new);

#if 0
    std::cout << l1 << " " << l2 << " " << l3 << std::endl;
    std::cout << l1_ <<  " " << l2_ << " " << l3_ << std::endl;
    std::cout << cA12 << " " << cA31 << " " << cA23 << " " << cA14 << " " << cA42 << std::endl;
    std::cout << cA12_ctrl << " " << cA31_ctrl << " " << cA23_ctrl << " " << cA14_ctrl << " " << cA42_ctrl << std::endl;
    std::cout << M[0] << std::endl;
    std::cout << coord << std::endl;
    std::cout << Point::GetDistance(M[0], coord) << std::endl;
#endif

    return Point::GetDistance(M[0], coord);
}

double test_solutions(unsigned iters, double lambda)
{
    Point M[5];
    double l, s;

    read_data(M, l, s);

    double d12 = Point::GetSquaredDistance(M[1], M[2]);
    double d23 = Point::GetSquaredDistance(M[2], M[3]);
    double d31 = Point::GetSquaredDistance(M[3], M[1]);

    double radius = get_radius(d12, d23, d31);

    std::vector<double> errors_coord(iters);

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0, lambda * radius / 3);


    for (auto & elem : errors_coord)
    {
        M[4] = M[1] + (M[2] - M[1]) * 0.5 + Point(distribution(generator),
                                                  distribution(generator),
                                                  distribution(generator));
        elem = solve(M, l, s);
    }

    double err_avg = std::accumulate(errors_coord.begin(), errors_coord.end(), 0.0) / iters;
    std::cout << lambda * radius / 3 << " " << err_avg  << " " << err_avg / radius << std::endl;
}

double test_solutions_exact(unsigned iters, double deviation)
{
    Point M[5];
    double l, s;

    read_data(M, l, s);

    double d12 = Point::GetSquaredDistance(M[1], M[2]);
    double d23 = Point::GetSquaredDistance(M[2], M[3]);
    double d31 = Point::GetSquaredDistance(M[3], M[1]);

    double radius = get_radius(d12, d23, d31);

    std::vector<double> errors_coord(iters);

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0, deviation * radius / 3);

    Point dist[] = { Point(),
                     M[1] - M[0],
                     M[2] - M[0],
                     M[3] - M[0],
                     M[4] - M[0]};

    // sides to top
    double l1 = dist[1].GetSquaredMagnitude();
    double l2 = dist[2].GetSquaredMagnitude();
    double l3 = dist[3].GetSquaredMagnitude();

    for (auto & elem : errors_coord)
    {
        double l[] = { 0,
                        l1 + distribution(generator),
                        l2 + distribution(generator),
                        l3 + distribution(generator) } ;
        auto p = find_coordinates(M, l);
        elem = Point::GetDistance(M[0], p);
    }

    double err_avg = std::accumulate(errors_coord.begin(), errors_coord.end(), 0.0) / iters;
    std::cout << deviation * radius / 3 << " " << err_avg << std::endl;
}


int main()
{
    Point M[5];
    double lambda = 0, sigma = 0;
    read_data(M, lambda, sigma);
    solve(M, lambda, sigma);

//    for (auto i = 0; i <= 100; ++i)
//    {
//        test_solutions(1000, 0.001 * i);
//    }

    for (auto i = 0; i <= 100; ++i)
    {
        test_solutions_exact(1000, 0.001 * i);
    }

    return 0;
}