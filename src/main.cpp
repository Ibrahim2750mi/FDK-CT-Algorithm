#include <cstdint>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <opencv2/opencv.hpp>
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

std::vector<std::vector<std::vector<double>>>
loadProjectionsOpenCV(const GeometryParams& params,
                      const std::string& basePath,
                      const std::string& ext = ".png")
{
    int A = params.numAngles;
    int R = params.numDetectorRows;
    int C = params.numDetectorCols;

    std::vector<std::vector<std::vector<double>>> projections(
        A, std::vector<std::vector<double>>(R, std::vector<double>(C, 0.0)));

    for (int a = 0; a < A; ++a) {
        std::ostringstream name;
        name << basePath << std::setw(3) << std::setfill('0') << a << ext;

        cv::Mat img = cv::imread(name.str(), cv::IMREAD_ANYDEPTH);
        if (img.empty()) {
            std::cerr << "ERROR: cannot open " << name.str() << std::endl;
            std::exit(1);
        }

        if (img.rows != R || img.cols != C) {
            std::cerr << "ERROR: image size mismatch in " << name.str() << std::endl;
            std::exit(1);
        }

        for (int r = 0; r < R; ++r)
            for (int c = 0; c < C; ++c)
                projections[a][r][c] = static_cast<double>(img.at<std::uint16_t>(r, c));
    }

    return projections;
}


int main() {
    // Parameters from Table 1 (Fig. 5 setup)
    GeometryParams params;
    params.d = 60.0;
    params.D = 120.0;
    params.numDetectorRows = 39;
    params.numDetectorCols = 65;
    params.detectorSpacing = 1;
    params.numAngles = 128;

    std::cout << "Creating phantom..." << std::endl;
    Phantom phantom;

    std::cout << "Creating reconstructor..." << std::endl;
    FDKReconstructor fdk(params);

    std::cout << "Generating projections..." << std::endl;
    auto projections = loadProjectionsOpenCV(params, );

    std::cout << "Reconstructing..." << std::endl;
    auto reconstruction = fdk.reconstruct(projections);

    std::cout << "Saving results..." << std::endl;
    // Save middle slice (z=0 plane)
    saveSlice(reconstruction[24], "midplane.txt");

    std::cout << "Done!" << std::endl;

    return 0;
}