#include "../include/fdk.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>

// ============================================================================
// EXACT Python Phantom Implementation
// ============================================================================
class PythonPhantom : public Phantom {
public:
    PythonPhantom() {
        std::cout << "PythonPhantom created (exact match to Python)" << std::endl;
    }

    double getDensity(double x, double y, double z) const override {
        // x, y, z are in NORMALIZED coordinates [-1, 1]
        double density = 0.0;

        // Central sphere (radius 0.5)
        double r1 = sqrt(x*x + y*y + z*z);
        if (r1 < 0.5) {
            density += 1.0 * (1.0 - (r1/0.5)*(r1/0.5));
        }

        // Right sphere (radius 0.2)
        double r2 = sqrt((x-0.6)*(x-0.6) + y*y + z*z);
        if (r2 < 0.2) {
            density += 0.8 * (1.0 - (r2/0.2)*(r2/0.2));
        }

        // Left sphere (radius 0.2)
        double r3 = sqrt((x+0.6)*(x+0.6) + y*y + z*z);
        if (r3 < 0.2) {
            density += 0.6 * (1.0 - (r3/0.2)*(r3/0.2));
        }

        // Top sphere (radius 0.15)
        double r4 = sqrt(x*x + (y-0.6)*(y-0.6) + z*z);
        if (r4 < 0.15) {
            density += 0.4 * (1.0 - (r4/0.15)*(r4/0.15));
        }

        // Bottom sphere (radius 0.15)
        double r5 = sqrt(x*x + (y+0.6)*(y+0.6) + z*z);
        if (r5 < 0.15) {
            density += 0.3 * (1.0 - (r5/0.15)*(r5/0.15));
        }

        return density;
    }
};

// ============================================================================
// Save results for Python comparison
// ============================================================================
void saveSlice(const std::vector<std::vector<double>>& slice,
               const std::string& filename) {
    std::ofstream file(filename);
    int N = slice.size();

    for (int y = 0; y < N; y++) {
        for (int x = 0; x < N; x++) {
            file << slice[y][x] << " ";
        }
        file << "\n";
    }
    file.close();
    std::cout << "Saved " << filename << " (" << N << "x" << N << ")" << std::endl;
}

// ============================================================================
// Debug: Test single ray through phantom
// ============================================================================
void testSingleRay(const PythonPhantom& phantom) {
    std::cout << "\n=== TESTING SINGLE RAY ===" << std::endl;

    // Simple ray through center
    double sx = 400.0, sy = 0.0, sz = 0.0;  // Source at (d, 0, 0)
    double dx = -400.0, dy = 0.0, dz = 0.0; // Detector at (-d, 0, 0)

    std::cout << "Source: (" << sx << ", " << sy << ", " << sz << ")" << std::endl;
    std::cout << "Detector: (" << dx << ", " << dy << ", " << dz << ")" << std::endl;

    // Ray direction
    double vx = dx - sx;
    double vy = dy - sy;
    double vz = dz - sz;
    double L = sqrt(vx*vx + vy*vy + vz*vz);
    vx /= L; vy /= L; vz /= L;

    std::cout << "Ray direction: (" << vx << ", " << vy << ", " << vz << ")" << std::endl;
    std::cout << "Ray length: " << L << " mm" << std::endl;

    // Sample along ray
    double sum = 0.0;
    int steps = 800;
    double dt = L / steps;

    std::cout << "\nSampling " << steps << " points along ray:" << std::endl;

    int nonzero_count = 0;
    for (int i = 0; i < steps; i++) {
        double t = i * dt;
        double x = sx + t * vx;
        double y = sy + t * vy;
        double z = sz + t * vz;

        // Normalize to [-1, 1]
        double xn = x / 100.0;
        double yn = y / 100.0;
        double zn = z / 100.0;

        double density = phantom.getDensity(xn, yn, zn);

        if (density > 0.0) {
            nonzero_count++;
            if (nonzero_count <= 5) {  // Print first 5 non-zero samples
                std::cout << "  Step " << i << ": pos=(" << x << "," << y << "," << z
                          << ") normalized=(" << xn << "," << yn << "," << zn
                          << ") density=" << density << std::endl;
            }
        }

        sum += density * dt;
    }

    std::cout << "Non-zero samples: " << nonzero_count << "/" << steps << std::endl;
    std::cout << "Line integral: " << sum << std::endl;
}

int main() {
    std::cout << "==========================================" << std::endl;
    std::cout << "DEBUG FDK TEST" << std::endl;
    std::cout << "==========================================" << std::endl;

    // Use small parameters for quick test
    GeometryParams params;
    params.d = 400.0;           // source → iso
    params.sdd = 800.0;         // source → detector

    double pitch_per_rotation = 40.0;   // mm
    params.h = pitch_per_rotation / (2 * M_PI);

    params.numDetectorCols = 256;
    params.numDetectorRows = 64;
    params.detectorSpacing = 2.0;
    params.numAngles = 90;
    params.reconSize = 64;
    params.reconSpacing = 3.125;

    params.helicalPitch = 2 * M_PI * params.h;
    params.detectorColSize = params.numDetectorCols * params.detectorSpacing;
    params.detectorRowSize = params.numDetectorRows * params.detectorSpacing;

    std::cout << "\nGeometry Parameters:" << std::endl;
    std::cout << "  d (SID): " << params.d << " mm" << std::endl;
    std::cout << "  sdd (SDD): " << params.sdd << " mm" << std::endl;
    std::cout << "  h: " << params.h << " mm/rad" << std::endl;
    std::cout << "  Detector: " << params.numDetectorCols << " x " << params.numDetectorRows << std::endl;
    std::cout << "  Detector spacing: " << params.detectorSpacing << " mm" << std::endl;
    std::cout << "  Recon volume: " << params.reconSize << "^3" << std::endl;
    std::cout << "  Recon spacing: " << params.reconSpacing << " mm" << std::endl;
    std::cout << "  FOV: " << params.reconSize * params.reconSpacing << " mm" << std::endl;

    // Create phantom
    PythonPhantom phantom;

    // Test phantom sampling
    std::cout << "\n=== TESTING PHANTOM SAMPLING ===" << std::endl;
    std::cout << "Center (0,0,0): " << phantom.getDensity(0,0,0) << std::endl;
    std::cout << "Right sphere center (0.6,0,0): " << phantom.getDensity(0.6,0,0) << std::endl;
    std::cout << "Left sphere center (-0.6,0,0): " << phantom.getDensity(-0.6,0,0) << std::endl;
    std::cout << "Outside (1,0,0): " << phantom.getDensity(1,0,0) << std::endl;

    // Test single ray
    testSingleRay(phantom);

    // Create reconstructor
    FDKReconstructor reconstructor(params);

    auto start = std::chrono::high_resolution_clock::now();

    // Generate projections
    std::cout << "\n=== GENERATING PROJECTIONS ===" << std::endl;
    auto projections = reconstructor.generateProjections(phantom);

    auto mid = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> proj_time = mid - start;
    std::cout << "Projection time: " << proj_time.count() << "s" << std::endl;

    // Detailed projection statistics
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
                if (val > 1e-10) nonzero_pixels++;
                total_pixels++;
            }
        }
    }

    std::cout << "Projection min: " << minProj << std::endl;
    std::cout << "Projection max: " << maxProj << std::endl;
    std::cout << "Projection mean: " << sumProj / total_pixels << std::endl;
    std::cout << "Non-zero pixels: " << nonzero_pixels << "/" << total_pixels
              << " (" << (100.0 * nonzero_pixels / total_pixels) << "%)" << std::endl;

    // Save first projection
    if (!projections.empty()) {
        std::ofstream projFile("first_projection.txt");
        int rows = projections[0].size();
        int cols = projections[0][0].size();

        std::cout << "\nFirst projection sample (center region):" << std::endl;
        int centerRow = rows / 2;
        int centerCol = cols / 2;
        for (int r = centerRow - 2; r <= centerRow + 2; r++) {
            for (int c = centerCol - 2; c <= centerCol + 2; c++) {
                std::cout << std::setw(10) << std::setprecision(4) << projections[0][r][c] << " ";
            }
            std::cout << std::endl;
        }

        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                projFile << projections[0][r][c] << " ";
            }
            projFile << "\n";
        }
        projFile.close();
        std::cout << "\nSaved first_projection.txt (" << rows << "x" << cols << ")" << std::endl;
    }

    if (maxProj < 1e-10) {
        std::cout << "\n*** ERROR: All projections are zero! ***" << std::endl;
        std::cout << "This indicates a problem with the forward projection geometry." << std::endl;
        return 1;
    }

    // Reconstruct
    std::cout << "\n=== RECONSTRUCTING ===" << std::endl;
    auto volume = reconstructor.reconstruct(projections);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> recon_time = end - mid;
    std::chrono::duration<double> total_time = end - start;

    std::cout << "\nTiming:" << std::endl;
    std::cout << "  Projections: " << proj_time.count() << "s" << std::endl;
    std::cout << "  Reconstruction: " << recon_time.count() << "s" << std::endl;
    std::cout << "  Total: " << total_time.count() << "s" << std::endl;

    // Volume statistics
    std::cout << "\n=== VOLUME STATISTICS ===" << std::endl;
    int N = volume.size();
    double minVal = 1e100, maxVal = -1e100, sum = 0.0;
    int nonzero_voxels = 0;

    for (const auto& slice : volume) {
        for (const auto& row : slice) {
            for (double val : row) {
                if (val < minVal) minVal = val;
                if (val > maxVal) maxVal = val;
                sum += val;
                if (fabs(val) > 1e-10) nonzero_voxels++;
            }
        }
    }

    std::cout << "Volume min: " << minVal << std::endl;
    std::cout << "Volume max: " << maxVal << std::endl;
    std::cout << "Volume mean: " << sum/(N*N*N) << std::endl;
    std::cout << "Non-zero voxels: " << nonzero_voxels << "/" << (N*N*N) << std::endl;

    // Sample center slice
    int midZ = N / 2;
    std::cout << "\nCenter XY slice sample:" << std::endl;
    for (int y = N/2 - 2; y <= N/2 + 2; y++) {
        for (int x = N/2 - 2; x <= N/2 + 2; x++) {
            std::cout << std::setw(10) << std::setprecision(4) << volume[midZ][y][x] << " ";
        }
        std::cout << std::endl;
    }

    // Save slices
    int midY = N / 2;
    int midX = N / 2;

    std::vector<std::vector<double>> xySlice(N, std::vector<double>(N));
    std::vector<std::vector<double>> xzSlice(N, std::vector<double>(N));
    std::vector<std::vector<double>> yzSlice(N, std::vector<double>(N));

    for (int y = 0; y < N; y++) {
        for (int x = 0; x < N; x++) {
            xySlice[y][x] = volume[midZ][y][x];
        }
    }

    for (int z = 0; z < N; z++) {
        for (int x = 0; x < N; x++) {
            xzSlice[z][x] = volume[z][midY][x];
        }
    }

    for (int z = 0; z < N; z++) {
        for (int y = 0; y < N; y++) {
            yzSlice[z][y] = volume[z][y][midX];
        }
    }

    saveSlice(xySlice, "cpp_xy_slice.txt");
    saveSlice(xzSlice, "cpp_xz_slice.txt");
    saveSlice(yzSlice, "cpp_yz_slice.txt");

    std::cout << "\nDone! Check the debug output above for issues." << std::endl;

    return 0;
}