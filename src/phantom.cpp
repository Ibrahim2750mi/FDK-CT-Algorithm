#include "../include/phantom.h"

Phantom::Phantom() {
    // From Table 1 in the paper
    ellipsoids = {
        {0, 0, 0, 40, 40, 1e9, 1.0},           // Outer cylinder
        {0, 0, 0, 34, 34, 1e9, -1.21},         // Inner cylinder (negative = less dense)
        {0, 0, 0, 30, 20, 20.98, 0.21},      // Central ellipsoid
        {-5, 0, 5, 10.95, 10.95, 10.95, 0.053},   // Small sphere
        {-7, -6, -5, 14.14, 16.73, 10.95, 0.316}, // Ellipsoid 5
        {8, 8, 2, 12, 8, 16, 0.158}          // Ellipsoid 6
    };
}

double Phantom::getDensity(double x, double y, double z) const {
    double total = 0.0;

    for (const auto& e : ellipsoids) {
        // Check if point is inside ellipsoid
        // Semi-axes are half the diameters
        double rx = e.dx / 2.0;
        double ry = e.dy / 2.0;
        double rz = e.dz / 2.0;

        double term = pow((x - e.x) / rx, 2) +
                      pow((y - e.y) / ry, 2) +
                      pow((z - e.z) / rz, 2);

        if (term <= 1.0) {
            total += e.density;
        }
    }

    return total;
}