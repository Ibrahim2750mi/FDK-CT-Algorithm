#include "../include/fdk.h"

#include <opencv2/opencv.hpp>

#include <cmath>
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <vector>
#include <string>
#include <stdexcept>

// ============================================================================
// Load projections from PNG images
// pattern: printf-style pattern, e.g. "projections/proj_%03d.png"
// numAngles: number of projection images / views
// Returns: projections[angle][row][col] as double
// ============================================================================
std::vector<std::vector<std::vector<double>>>
loadProjectionsFromPNGs(const std::string& pattern,
                        int numAngles,
                        int& outRows,
                        int& outCols)
{
    std::vector<std::vector<std::vector<double>>> projections;
    projections.reserve(numAngles);

    outRows = 0;
    outCols = 0;

    char filename[512];

    for (int a = 0; a < numAngles; ++a) {
        std::snprintf(filename, sizeof(filename), pattern.c_str(), a);
        std::string fname(filename);

        std::cout << "Loading projection " << a << ": " << fname << std::endl;

        // Read as grayscale (keep depth)
        cv::Mat img = cv::imread(fname, cv::IMREAD_UNCHANGED);
        if (img.empty()) {
            throw std::runtime_error("Failed to load image: " + fname);
        }

        if (img.channels() > 1) {
            cv::cvtColor(img, img, cv::COLOR_BGR2GRAY);
        }

        // Convert to double
        cv::Mat img64;
        img.convertTo(img64, CV_64F);

        if (a == 0) {
            outRows = img64.rows;
            outCols = img64.cols;
            std::cout << "Detected detector size from first image: "
                      << outCols << " x " << outRows << std::endl;
        } else {
            if (img64.rows != outRows || img64.cols != outCols) {
                throw std::runtime_error("Inconsistent image size at angle "
                                         + std::to_string(a));
            }
        }

        // Copy into std::vector structure [row][col]
        std::vector<std::vector<double>> angle(outRows,
                                               std::vector<double>(outCols, 0.0));
        for (int r = 0; r < outRows; ++r) {
            const double* rowPtr = img64.ptr<double>(r);
            for (int c = 0; c < outCols; ++c) {
                angle[r][c] = rowPtr[c];
            }
        }

        projections.push_back(std::move(angle));
    }

    std::cout << "Loaded " << projections.size() << " projections." << std::endl;
    return projections;
}

// ============================================================================
// Save full 3D volume as RAW (float32), layout: z-major, then y, then x
// volume[z][y][x]
// ============================================================================
void saveRawVolume(const std::vector<std::vector<std::vector<double>>>& volume,
                   const std::string& filename)
{
    if (volume.empty() || volume[0].empty() || volume[0][0].empty()) {
        std::cerr << "Volume is empty, not saving RAW." << std::endl;
        return;
    }

    int Nz = static_cast<int>(volume.size());
    int Ny = static_cast<int>(volume[0].size());
    int Nx = static_cast<int>(volume[0][0].size());

    std::ofstream out(filename, std::ios::binary);
    if (!out) {
        std::cerr << "Failed to open " << filename << " for writing." << std::endl;
        return;
    }

    for (int z = 0; z < Nz; ++z) {
        for (int y = 0; y < Ny; ++y) {
            for (int x = 0; x < Nx; ++x) {
                float v = static_cast<float>(volume[z][y][x]);
                out.write(reinterpret_cast<const char*>(&v), sizeof(float));
            }
        }
    }

    out.close();
    std::cout << "Saved RAW volume to " << filename
              << " (" << Nz << " x " << Ny << " x " << Nx
              << ", float32)" << std::endl;
}

// ============================================================================
// Save all axial slices (z = const planes) as PNG
// Filenames: pattern with %03d, e.g. "axial_%03d.png"
// Each slice is normalized independently to [0,255]
// ============================================================================
void saveAxialSlicesAsPNG(const std::vector<std::vector<std::vector<double>>>& volume,
                          const std::string& pattern)
{
    if (volume.empty() || volume[0].empty() || volume[0][0].empty()) {
        std::cerr << "Volume is empty, not saving PNG slices." << std::endl;
        return;
    }

    int Nz = static_cast<int>(volume.size());
    int Ny = static_cast<int>(volume[0].size());
    int Nx = static_cast<int>(volume[0][0].size());

    char filename[256];

    for (int z = 0; z < Nz; ++z) {
        cv::Mat img(Ny, Nx, CV_32F);

        // Copy slice z into img
        for (int y = 0; y < Ny; ++y) {
            float* rowPtr = img.ptr<float>(y);
            for (int x = 0; x < Nx; ++x) {
                rowPtr[x] = static_cast<float>(volume[z][y][x]);
            }
        }

        // Normalize to [0,255] for visualization
        double mn, mx;
        cv::minMaxLoc(img, &mn, &mx);

        cv::Mat imgNorm;
        if (mx > mn) {
            imgNorm = (img - mn) / (mx - mn);
            imgNorm *= 255.0;
        } else {
            imgNorm = cv::Mat::zeros(img.size(), CV_32F);
        }

        cv::Mat img8;
        imgNorm.convertTo(img8, CV_8U);

        std::snprintf(filename, sizeof(filename), pattern.c_str(), z);
        std::string fname(filename);

        if (!cv::imwrite(fname, img8)) {
            std::cerr << "Failed to write PNG slice: " << fname << std::endl;
        } else {
            std::cout << "Saved axial slice " << z << " -> " << fname << std::endl;
        }
    }

    std::cout << "Saved " << Nz << " axial PNG slices." << std::endl;
}

int main(int argc, char** argv) {
    std::cout << "==========================================" << std::endl;
    std::cout << "HELICAL FDK RECONSTRUCTION (PNG INPUT)" << std::endl;
    std::cout << "==========================================" << std::endl;

    if (argc < 3) {
        std::cout << "Usage: " << argv[0]
                  << " <projection_pattern> <num_angles>\n"
                  << "Example: " << argv[0]
                  << " projections/proj_%03d.png 90\n";
        return 1;
    }

    std::string projectionPattern = argv[1];
    int numAngles = std::stoi(argv[2]);

    // -------------------------------
    // Set geometry parameters
    // -------------------------------
    GeometryParams params;

    // Basic helical geometry (you can tune these to match your scanner / sim)
    params.d   = 400.0;  // source-to-iso (SID)
    params.sdd = 800.0;  // source-to-detector (SDD)

    double pitch_per_rotation = 40.0;   // mm per 2π
    params.h = pitch_per_rotation / (2 * M_PI);  // mm / rad

    params.numAngles = numAngles;

    // Detector spacing (mm) – adjust to match your data
    params.detectorSpacing = 2.0;

    // Reconstruction volume
    params.reconSize    = 64;
    params.reconSpacing = 3.125;

    // We'll set numDetectorRows/Cols after we inspect first image.
    int detRows = 0, detCols = 0;

    // -------------------------------
    // Load projections from PNGs
    // -------------------------------
    auto io_start = std::chrono::high_resolution_clock::now();

    std::vector<std::vector<std::vector<double>>> projections;
    try {
        projections = loadProjectionsFromPNGs(projectionPattern,
                                              numAngles,
                                              detRows,
                                              detCols);
    } catch (const std::exception& e) {
        std::cerr << "Error while loading projections: " << e.what() << std::endl;
        return 1;
    }

    auto io_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> io_time = io_end - io_start;
    std::cout << "Projection loading time: " << io_time.count() << " s" << std::endl;

    // Now we know detector size
    params.numDetectorRows = detRows;
    params.numDetectorCols = detCols;

    params.helicalPitch    = 2 * M_PI * params.h;
    params.detectorColSize = params.numDetectorCols * params.detectorSpacing;
    params.detectorRowSize = params.numDetectorRows * params.detectorSpacing;

    std::cout << "\nGeometry Parameters:" << std::endl;
    std::cout << "  d (SID): " << params.d << " mm" << std::endl;
    std::cout << "  sdd (SDD): " << params.sdd << " mm" << std::endl;
    std::cout << "  h: " << params.h << " mm/rad" << std::endl;
    std::cout << "  Helical pitch (per rotation): " << params.helicalPitch << " mm" << std::endl;
    std::cout << "  Detector: " << params.numDetectorCols << " x " << params.numDetectorRows << std::endl;
    std::cout << "  Detector spacing: " << params.detectorSpacing << " mm" << std::endl;
    std::cout << "  Recon volume: " << params.reconSize << "^3" << std::endl;
    std::cout << "  Recon spacing: " << params.reconSpacing << " mm" << std::endl;
    std::cout << "  FOV: " << params.reconSize * params.reconSpacing << " mm" << std::endl;

    // -------------------------------
    // Projection statistics
    // -------------------------------
    std::cout << "\n=== PROJECTION STATISTICS ===" << std::endl;
    double minProj = 1e100, maxProj = -1e100, sumProj = 0.0;
    int nonzero_pixels = 0;
    int total_pixels = 0;

    for (const auto& angle : projections) {
        for (const auto& row : angle) {
            for (double val : row) {
                minProj = std::min(minProj, val);
                maxProj = std::max(maxProj, val);
                sumProj += val;
                if (std::fabs(val) > 1e-10) nonzero_pixels++;
                total_pixels++;
            }
        }
    }

    std::cout << "Projection min: " << minProj << std::endl;
    std::cout << "Projection max: " << maxProj << std::endl;
    std::cout << "Projection mean: " << sumProj / total_pixels << std::endl;
    std::cout << "Non-zero pixels: " << nonzero_pixels << "/" << total_pixels
              << " (" << (100.0 * nonzero_pixels / total_pixels) << "%)" << std::endl;

    if (maxProj < 1e-10) {
        std::cout << "\n*** ERROR: All projections are ~zero! ***" << std::endl;
        std::cout << "Check that the PNG images contain the expected data "
                  << "(log-intensities / line integrals)." << std::endl;
        return 1;
    }

    // -------------------------------
    // Reconstruction
    // -------------------------------
    std::cout << "\n=== RECONSTRUCTING ===" << std::endl;
    auto recon_start = std::chrono::high_resolution_clock::now();

    FDKReconstructor reconstructor(params);
    auto volume = reconstructor.reconstruct(projections);

    auto recon_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> recon_time = recon_end - recon_start;

    std::cout << "\nTiming:" << std::endl;
    std::cout << "  Projection loading: " << io_time.count() << " s" << std::endl;
    std::cout << "  Reconstruction:     " << recon_time.count() << " s" << std::endl;
    std::cout << "  Total:              " << (io_time + recon_time).count() << " s" << std::endl;

    // -------------------------------
    // Volume statistics
    // -------------------------------
    std::cout << "\n=== VOLUME STATISTICS ===" << std::endl;
    int N = static_cast<int>(volume.size());
    double minVal = 1e100, maxVal = -1e100, sum = 0.0;
    int nonzero_voxels = 0;

    for (const auto& slice : volume) {
        for (const auto& row : slice) {
            for (double val : row) {
                if (val < minVal) minVal = val;
                if (val > maxVal) maxVal = val;
                sum += val;
                if (std::fabs(val) > 1e-10) nonzero_voxels++;
            }
        }
    }

    std::cout << "Volume min: " << minVal << std::endl;
    std::cout << "Volume max: " << maxVal << std::endl;
    std::cout << "Volume mean: " << sum / (N * N * N) << std::endl;
    std::cout << "Non-zero voxels: " << nonzero_voxels
              << "/" << (N * N * N) << std::endl;

    // -------------------------------
    // Save outputs: RAW volume + axial PNG slices
    // -------------------------------
    saveRawVolume(volume, "recon_volume.raw");
    saveAxialSlicesAsPNG(volume, "axial_%03d.png");

    std::cout << "\nDone! Outputs:\n"
              << "  - recon_volume.raw (float32, z-y-x)\n"
              << "  - axial_###.png (axial slices)\n";
    return 0;
}
