#ifndef PHANTOM_H
#define PHANTOM_H

#include <vector>
#include <cmath>

struct Ellipsoid {
    double x, y, z;           // Center position
    double dx, dy, dz;        // Diameters (semi-axes will be half of these)
    double density;
};

class Phantom {
public:
    Phantom();
    virtual Phantom(double x, double y, double z) const;

private:
    std::vector<Ellipsoid> ellipsoids;
};

#endif