#ifndef FDK_H
#define FDK_H

#include <math.h>
#include <vector>
#include <string>

struct GeometryParams {
    // Source-detector geometry
    double d = 400.0;          // Source-to-isocenter distance (mm)
    double sdd = 800.0;        // Source-to-detector distance (mm)
    double h = 10.0;           // Helical pitch parameter (mm/rad)
    double helicalPitch;       // 2πh

    // Detector parameters
    int numDetectorCols = 512;
    int numDetectorRows = 64;
    double detectorSpacing = 1.0;  // mm
    double detectorColSize;         // cols * spacing
    double detectorRowSize;         // rows * spacing

    // Acquisition parameters
    int numAngles = 360;           // Number of projection angles

    // Reconstruction parameters
    int reconSize = 256;           // Reconstruction volume size (cubic)
    double reconSpacing = 1.0;     // mm

    GeometryParams() {
        helicalPitch = 2 * M_PI * h;
        detectorColSize = numDetectorCols * detectorSpacing;
        detectorRowSize = numDetectorRows * detectorSpacing;
    }
};

class Phantom {
public:
    virtual double getDensity(double x, double y, double z) const = 0;
    virtual ~Phantom() {}
};

class FDKReconstructor {
public:
    FDKReconstructor(const GeometryParams& params);

    // Generate helical cone-beam projections from phantom
    std::vector<std::vector<std::vector<double>>>
    generateProjections(const Phantom& phantom);

    // Main reconstruction function
    std::vector<std::vector<std::vector<double>>>
    reconstruct(const std::vector<std::vector<std::vector<double>>>& projections);

private:
    GeometryParams params_;

    // Helper functions
    void cosineWeight(std::vector<std::vector<double>>& projection);
    void applyParkerWeight(std::vector<std::vector<std::vector<double>>>& projections);
    void filterProjection(std::vector<std::vector<double>>& projection);
    void backproject(const std::vector<std::vector<std::vector<double>>>& filteredProjections,
                     std::vector<std::vector<std::vector<double>>>& volume) const;

    double computeLineIntegral(const Phantom& phantom,
                               double srcX, double srcY, double srcZ,
                               double detX, double detY, double detZ) const;

    // FDK-2 specific functions
    double getTiltAngle() const { return atan2(params_.h, params_.d); }

    // Parker window helper
    double parkerWindow(double beta, double gamma) const;

    // Filter helper
    std::vector<double> rampFilter(const std::vector<double>& signal) const;
};

#endif