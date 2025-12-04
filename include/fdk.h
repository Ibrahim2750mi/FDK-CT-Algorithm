#ifndef FDK_H
#define FDK_H

#include <vector>
#include <cmath>
#include "phantom.h"

struct GeometryParams {
    double d;           // Source to rotation axis distance
    double D;           // Source to detector distance
    int numDetectorRows;    // Number of detector rows (Z direction)
    int numDetectorCols;    // Number of detector columns (Y direction)
    double detectorSpacing; // Spacing between detector pixels
    int numAngles;          // Number of projection angles

    double h;
};

class FDKReconstructor {
public:
    FDKReconstructor(const GeometryParams& params);

    // Generate projection data from phantom
    std::vector<std::vector<std::vector<double>>> generateProjections(
        const Phantom& phantom);

    // Main reconstruction function
    std::vector<std::vector<std::vector<double>>> reconstruct(
        const std::vector<std::vector<std::vector<double>>>& projections);

private:
    GeometryParams params_;

    // Cosine weighting (Section 3, before Eq. 31)
    void cosineWeight(std::vector<std::vector<double>>& projection);

    // Filtering (Eq. 31-33)
    void filterProjection(std::vector<std::vector<double>>& projection);

    // Backprojection (Eq. 28)
    void backproject(
        const std::vector<std::vector<std::vector<double>>>& filteredProjections,
        std::vector<std::vector<std::vector<double>>>& volume) const;

    // compute line integral through phantom
    double computeLineIntegral(const Phantom& phantom,
                               double srcX, double srcY, double srcZ,
                               double detX, double detY, double detZ);
};

#endif