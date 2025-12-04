#include <iostream>
#include <fstream>
#include "../include/fdk.h"
#include "../include/phantom.h"

void saveSlice(const std::vector<std::vector<double>>& slice,
               const std::string& filename) {
    std::ofstream file(filename);
    for (const auto& row : slice) {
        for (double val : row) {
            file << val << " ";
        }
        file << "\n";
    }
}

int main() {
    // Parameters from Table 1 (Fig. 5 setup)
    GeometryParams params;
    params.d = 60.0;
    params.D = 60.0;
    params.numDetectorRows = 39;
    params.numDetectorCols = 65;
    params.detectorSpacing = 1;
    params.numAngles = 32;

    params.h = 46.28 / (2.0 * M_PI);   // Pitch P = 46.28mm, so h = P/(2π)

    std::cout << "Helical pitch P = " << (2 * M_PI * params.h) << std::endl;

    std::cout << "Creating phantom..." << std::endl;
    Phantom phantom;

    std::cout << "Creating reconstructor..." << std::endl;
    FDKReconstructor fdk(params);

    std::cout << "Generating helical projections..." << std::endl;
    auto projections = fdk.generateProjections(phantom);

    // Debug: check projection values
    double maxP = -1e10, minP = 1e10;
    for (const auto& proj : projections) {
        for (const auto& row : proj) {
            for (double val : row) {
                maxP = std::max(maxP, val);
                minP = std::min(minP, val);
            }
        }
    }
    std::cout << "Projection range: [" << minP << ", " << maxP << "]" << std::endl;

    std::cout << "Reconstructing..." << std::endl;
    auto reconstruction = fdk.reconstruct(projections);

    std::cout << "Saving results..." << std::endl;
    saveSlice(reconstruction[24], "midplane.txt");

    // Save multiple slices to see helical effect
    saveSlice(reconstruction[0], "slice_bottom.txt");
    saveSlice(reconstruction[12], "slice_lower.txt");
    saveSlice(reconstruction[24], "slice_mid.txt");
    saveSlice(reconstruction[36], "slice_upper.txt");
    saveSlice(reconstruction[48], "slice_top.txt");

    std::cout << "Done!" << std::endl;
    return 0;
}