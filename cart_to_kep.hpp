#ifndef MY_PROJECT_CART_TO_KEP_HPP
#define MY_PROJECT_CART_TO_KEP_HPP
// Cartesian to keplerian
#include <iostream>
#include <array>
#include <cmath>
using namespace std;
const double M_Eath = 5.972e24; //Kg

std::array<double,6> cart_to_kep(const std::array<double,6>& temp) {
    double t = 0;
    //  x, y, z, x', y', z'
    std::array<double, 6> cart_arr = {0, 0, 0, 0, 0, 0};
    for (int i = 0; i < 6; ++i) {
        cart_arr[i] = temp[i];
    }

    double x = cart_arr[0];
    double y = cart_arr[1];
    double z = cart_arr[2];
    double x1 = cart_arr[3];
    double y1 = cart_arr[4];
    double z1 = cart_arr[5];

    // Computing the specific angular momentum and checking for a degenerate orbit but not yet
    double h1 = y * z1 - z * y1;
    double h2 = -x * z1 + z * x1;
    double h3 = x * y1 - x1 * y;
    double h = sqrt(std::pow(h1, 2) + std::pow(h2, 2) + std::pow(h3, 2));

    // Computing the radius, r, and velocity, v
    double r = sqrt(std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2));
    double v = sqrt(std::pow(x1, 2) + std::pow(y1, 2) + std::pow(z1, 2));

    // Computing the specific energy, E, and verify elliptical motion
    double mu = 6.6743 * 10e-11 * M_Eath;
    double E = std::pow(v, 2) / 2 - mu / r;

    //Computing semi-major axis, a
    double a = -mu / (2 * E);

    // Computing eccentricity, e
    double e = sqrt(1 - std::pow(h, 2) / (a * mu));

    // Computing inclination, i, (0° → 180°)
    double i = std::acos(h3 / h);

    // Compute right ascension of the ascending node, Ω, (0° → 360°)
    double omega = std::atan2(h1, -h2);

    // Compute argument of latitude, ω +ν , (0° → 360°)
    double l_om_plus_nu = atan2(z / sin(i), (x * cos(omega) + y * sin(omega)));

    // Computing true anomaly, ν , (0° → 360°)
    double p = a * (1 - std::pow(e, 2));
    double nu = atan2(std::pow(p / mu, 1 / 2) * (x1 * x + y1 * y + z1 * z), p - r);

    // Computing argument of periapse, ω, (0° → 360°)
    double l_om = l_om_plus_nu - nu;

    // Computing eccentric anomaly, EA, (0° → 360°)
    double EA = 2 * atan(std::pow((1 - e) / (1 + e), 1 / 2) * tan(nu / 2));

    // Computing the time of periapse passage, T, EA is acceptable only in radians
    double n = std::pow(mu / std::pow(a, 3), 1 / 2);
    double b_T = t - 1 * (EA - e * sin(EA)) / n;//
    std::array<double,6> result_in_kep{a,e,i, l_om, omega, EA};
    return  result_in_kep;
}
#endif //MY_PROJECT_CART_TO_KEP_HPP
