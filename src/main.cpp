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
    params.detectorSpacing = 1.0;
    params.numAngles = 32;

    std::cout << "Creating phantom..." << std::endl;
    Phantom phantom;

    std::cout << "Creating reconstructor..." << std::endl;
    FDKReconstructor fdk(params);

    std::cout << "Generating projections..." << std::endl;
    auto projections = fdk.generateProjections(phantom);

    std::cout << "Reconstructing..." << std::endl;
    auto reconstruction = fdk.reconstruct(projections);

    std::cout << "Saving results..." << std::endl;
    // Save middle slice (z=0 plane)
    saveSlice(reconstruction[24], "midplane.txt");

    std::cout << "Done!" << std::endl;

    return 0;
}