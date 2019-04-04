#include <iostream>

#define M_PI           3.14159265358979323846  /* pi */
#define DEG_TO_RAD(x)  (180 / M_PI * (x))

#include "point.hpp"

enum Ind
{
    I1_2 = 0,
    I1_3,
    I2_3
};

int inverse_third_order_matrix(double matrix[3][3])
{
    double a = matrix[0][0], b = matrix[0][1], c = matrix[0][2];
    double d = matrix[1][0], e = matrix[1][1], f = matrix[1][2];
    double g = matrix[2][0], h = matrix[2][1], i = matrix[2][2];

    double A = (e * i - f * h), B = -(d * i - f * g);
    double C = (d * h - e * g), D = -(b * i - c * h);
    double E = (a * i - c * g), F = -(a * h - b * g);
    double G = (b * f - c * e), H = -(a * f - c * d);
    double I = (a * e - b * d);

    double det = a * A + b * B + c * C;

    if (det == 0)
    {
        std::cerr << "Matrix is singular" << std::endl;
        return EXIT_FAILURE;
    }

    matrix[0][0] = A / det; matrix[0][1] = D / det; matrix[0][2] = G / det;
    matrix[1][0] = B / det; matrix[1][1] = E / det; matrix[1][2] = H / det;
    matrix[2][0] = C / det; matrix[2][1] = F / det; matrix[2][2] = I / det;

    return EXIT_SUCCESS;
}

int get_df_1(const double l[4], const double cos[3], double DF[3][3])
{
    double a = 2 * (l[1] - l[2] * cos[I1_2]), b = 2 * (l[2] - l[1] * cos[I1_2]);
    double c = 2 * (l[1] - l[3] * cos[I1_3]), d = 2 * (l[3] - l[1] * cos[I1_3]);
    double e = 2 * (l[2] - l[3] * cos[I2_3]), f = 2 * (l[3] - l[3] * cos[I2_3]);

    double delta = -(a * d * e + b * c * f);
    if (0 == delta)
    {
        std::cerr << "DF is singular!" << std::endl;
        return EXIT_FAILURE;
    }

    DF[0][0] = -d * e / delta; DF[0][1] = -b * f / delta; DF[0][2] =  b * d / delta;
    DF[1][0] = -c * f / delta; DF[1][1] =  a * f / delta; DF[1][2] = -a * d / delta;
    DF[2][0] =  c * e / delta; DF[2][1] = -a * e / delta; DF[2][2] = -b * c / delta;

    return EXIT_SUCCESS;
}

void read_data(Point M[4], Point &M0_start, Point &rotation, double &eps, unsigned &nmax)
{
    // TODO: read data from the file
    double phi = 35, psi = 270, theta = 5;

    M[0] = Point( 0.5, -0.5, 2.0);
    M[1] = Point( 2.0,  1.0, 0.0);
    M[2] = Point(-2.0,  1.0, 0.0);
    M[3] = Point( 1.0, -2.0, 0.0);

    M0_start = M[0] + Point(0.2, -0.1, 0.3);

    rotation = Point(DEG_TO_RAD(phi), DEG_TO_RAD(psi), DEG_TO_RAD(theta));

    eps = 1e-10;
    nmax = 10;
}

int solve(const Point M[4], const Point &M0_start, Point rotation, double eps, unsigned n_max)
{
    Point dist[4] = { Point(),
                      M[1] - M[0],
                      M[2] - M[0],
                      M[3] - M[0] };
    auto M0_new = M0_start;

    double l[4] = { 0,
                    dist[1].GetMagnitude(),
                    dist[2].GetMagnitude(),
                    dist[3].GetMagnitude() };
    double l_new[4] = { 0,
                        (M[1] - M0_start).GetMagnitude(),
                        (M[2] - M0_start).GetMagnitude(),
                        (M[3] - M0_start).GetMagnitude()};

    double d[3] = { (M[2] - M[1]).GetMagnitude(),
                    (M[3] - M[1]).GetMagnitude(),
                    (M[3] - M[2]).GetMagnitude() };
    double cos_ctrl[3] = { dist[1].GetAngle(dist[2]),
                           dist[1].GetAngle(dist[3]),
                           dist[2].GetAngle(dist[3])};
    double cos[3] = { (l[1] * l[1] + l[2] * l[2] - d[I1_2] * d[I1_2]) / (2 * l[1] * l[2]),
                      (l[1] * l[1] + l[3] * l[3] - d[I1_3] * d[I1_3]) / (2 * l[1] * l[3]),
                      (l[2] * l[2] + l[3] * l[3] - d[I2_3] * d[I2_3]) / (2 * l[2] * l[3]) };
    double DF[3][3] = { 0 };

    Point dist_rotated[] = { M[0],
                             dist[1].GetRotated(rotation),
                             dist[2].GetRotated(rotation),
                             dist[3].GetRotated(rotation) };
    double l_rotated[] = { 0,
                           dist_rotated[1].GetMagnitude(),
                           dist_rotated[2].GetMagnitude(),
                           dist_rotated[3].GetMagnitude()};
    double cos_rotated[] = { dist_rotated[1].GetAngle(dist_rotated[2]),
                             dist_rotated[1].GetAngle(dist_rotated[3]),
                             dist_rotated[2].GetAngle(dist_rotated[3])};

    std::cout << "DIST:     " << dist[1] << " " << dist[2] << " " << dist[3] << std::endl;
    std::cout << "DIST_rot: " << dist_rotated[1] << " " << dist_rotated[2] << " " << dist_rotated[3] << std::endl;
    std::cout << "L:        " << l[1] << " " << l[2] << " " << l[3] << std::endl;
    std::cout << "L_rot:    " << l_rotated[1] << " " << l_rotated[2] << " " << l_rotated[3] << std::endl;
    std::cout << "Cos:      " << cos[I1_2] << " " << cos[I1_3] << " " << cos[I2_3] << std::endl;
    std::cout << "Cos_rot:  " << cos_rotated[I1_2] << " " << cos_rotated[I1_3] << " " << cos_rotated[I2_3] << std::endl;

    // Step 1: Find side lenghts by base point coordinates
    auto err = Point::GetDistance(Point(l[1], l[2], l[3]), Point(l_new[1], l_new[2], l_new[3]));
    unsigned iter = 0;
    for (iter = 0; iter < n_max && err > eps; ++iter)
    {
        double l_old[] = { 0, l_new[1], l_new[2], l_new[3]};

        double f1 = l_old[1] * l_old[1] + l_old[2] * l_old[2] - 2 * l_old[1] * l_old[2] * cos_ctrl[I1_2] - d[I1_2] * d[I1_2];
        double f2 = l_old[1] * l_old[1] + l_old[3] * l_old[3] - 2 * l_old[1] * l_old[3] * cos_ctrl[I1_3] - d[I1_3] * d[I1_3];
        double f3 = l_old[2] * l_old[2] + l_old[3] * l_old[3] - 2 * l_old[2] * l_old[3] * cos_ctrl[I2_3] - d[I2_3] * d[I2_3];

        get_df_1(l_old, cos_ctrl, DF);

        l_new[1] = l_old[1] - DF[0][0] * f1 - DF[0][1] * f2 - DF[0][2] * f3;
        l_new[2] = l_old[2] - DF[1][0] * f1 - DF[1][1] * f2 - DF[1][2] * f3;
        l_new[3] = l_old[3] - DF[2][0] * f1 - DF[2][1] * f2 - DF[2][2] * f3;

        err = Point::GetDistance(Point(l_old[1], l_old[2], l_old[3]), Point(l_new[1], l_new[2], l_new[3]));
    }

    std::cout << "error: " << err << " " << "iter: " << iter << std::endl;
    std::cout << "l:     " << l[1] << " " << l[2] << " " << l[3] << std::endl;
    std::cout << "l_new: " << l_new[1] << " " << l_new[2] << " " << l_new[3] << std::endl;

    // Step 2: Find coordinate of the apex based on side lengths
    err = Point::GetDistance(M0_start, M[0]);
    for (iter = 0; iter < n_max && err > eps; ++iter)
    {
        Point M0_old = M0_new;

        double f1 = Point::GetSquaredDistance(M0_old, M[1]) - l_new[1] * l_new[1];
        double f2 = Point::GetSquaredDistance(M0_old, M[2]) - l_new[2] * l_new[2];
        double f3 = Point::GetSquaredDistance(M0_old, M[3]) - l_new[3] * l_new[3];

        dist[1] = M0_old - M[1];
        dist[2] = M0_old - M[2];
        dist[3] = M0_old - M[3];

        DF[0][0] = 2 * dist[1].GetX(); DF[0][1] = 2 * dist[1].GetY(); DF[0][2] = 2 * dist[1].GetZ();
        DF[1][0] = 2 * dist[2].GetX(); DF[1][1] = 2 * dist[2].GetY(); DF[1][2] = 2 * dist[2].GetZ();
        DF[2][0] = 2 * dist[3].GetX(); DF[2][1] = 2 * dist[3].GetY(); DF[2][2] = 2 * dist[3].GetZ();

        inverse_third_order_matrix(DF);

        Point fn = Point(DF[0][0] * f1 + DF[0][1] * f2 + DF[0][2] * f3,
                         DF[1][0] * f1 + DF[1][1] * f2 + DF[1][2] * f3,
                         DF[2][0] * f1 + DF[2][1] * f2 + DF[2][2] * f3);
        M0_new = M0_old - fn;

        err = Point::GetDistance(M0_new, M0_old);
    }

    std::cout << "error:  " << err << " " << "iter: " << iter << std::endl;
    std::cout << "M0:     " << M[0] << std::endl;
    std::cout << "M0_new: " << M0_new << std::endl;

    return 0;
}

int main()
{
    Point     M[4];
    Point M0_start = {};
    Point rotation;

    double eps;
    unsigned n_max;

    read_data(M, M0_start, rotation, eps, n_max);
    solve(M, M0_start, rotation, eps, n_max);

    return 0;
}
