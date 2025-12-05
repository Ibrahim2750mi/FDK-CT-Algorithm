// FDK IMPLEMENTATION WITH EXTENSIVE DEBUG OUTPUT

#include "../include/fdk.h"
#include <cmath>
#include <algorithm>
#include <complex>
#include <fftw3.h>
#include <iostream>
#include <iomanip>

const double PI = 3.14159265358979323846;

FDKReconstructor::FDKReconstructor(const GeometryParams& params)
    : params_(params) {
    std::cout << "\n=== FDK Reconstructor Initialized ===" << std::endl;
    std::cout << "d (source-to-axis): " << params_.d << std::endl;
    std::cout << "D (source-to-detector): " << params_.D << std::endl;
    std::cout << "Detector rows x cols: " << params_.numDetectorRows << " x " << params_.numDetectorCols << std::endl;
    std::cout << "Detector spacing: " << params_.detectorSpacing << std::endl;
    std::cout << "Number of angles: " << params_.numAngles << std::endl;

    if (params_.D <= params_.d) {
        std::cout << "WARNING: D <= d! Detector is at or behind rotation axis!" << std::endl;
    }
    std::cout << "Detector radius from axis: " << (params_.D - params_.d) << std::endl;
}

std::vector<std::vector<std::vector<double>>>
FDKReconstructor::generateProjections(const Phantom& phantom) {
    std::cout << "\n=== Generating Projections ===" << std::endl;

    std::vector<std::vector<std::vector<double>>> projections(
        params_.numAngles,
        std::vector<std::vector<double>>(
            params_.numDetectorRows,
            std::vector<double>(params_.numDetectorCols, 0.0)
        )
    );

    double minProj = 1e10, maxProj = -1e10;
    double sumProj = 0.0;
    int countNonZero = 0;

    for (int angleIdx = 0; angleIdx < params_.numAngles; angleIdx++) {
        double phi = 2.0 * PI * angleIdx / params_.numAngles;
        double srcX = params_.d * cos(phi);
        double srcY = params_.d * sin(phi);
        double srcZ = 0.0;

        for (int row = 0; row < params_.numDetectorRows; row++) {
            for (int col = 0; col < params_.numDetectorCols; col++) {
                double Y = (col - params_.numDetectorCols/2.0) * params_.detectorSpacing;
                double Z = (row - params_.numDetectorRows/2.0) * params_.detectorSpacing;

                double detRadius = params_.D - params_.d;
                double detCenterX = -detRadius * cos(phi);
                double detCenterY = -detRadius * sin(phi);
                double detCenterZ = 0.0;

                double uX = -sin(phi);
                double uY = cos(phi);

                double detX = detCenterX + Y * uX;
                double detY = detCenterY + Y * uY;
                double detZ = detCenterZ + Z;

                double val = computeLineIntegral(phantom, srcX, srcY, srcZ, detX, detY, detZ);
                projections[angleIdx][row][col] = val;

                minProj = std::min(minProj, val);
                maxProj = std::max(maxProj, val);
                sumProj += val;
                if (val > 1e-6) countNonZero++;
            }
        }
    }

    std::cout << "Projection statistics:" << std::endl;
    std::cout << "  Min: " << minProj << std::endl;
    std::cout << "  Max: " << maxProj << std::endl;
    std::cout << "  Mean: " << sumProj / (params_.numAngles * params_.numDetectorRows * params_.numDetectorCols) << std::endl;
    std::cout << "  Non-zero pixels: " << countNonZero << " / "
              << (params_.numAngles * params_.numDetectorRows * params_.numDetectorCols) << std::endl;

    // Debug: print one projection profile
    std::cout << "\nCenter row of first projection (cols 25-40):" << std::endl;
    int centerRow = params_.numDetectorRows / 2;
    for (int col = 25; col < 40; col++) {
        std::cout << std::setw(10) << std::fixed << std::setprecision(4)
                  << projections[0][centerRow][col] << " ";
    }
    std::cout << std::endl;

    return projections;
}

double FDKReconstructor::computeLineIntegral(
    const Phantom& phantom,
    double srcX, double srcY, double srcZ,
    double detX, double detY, double detZ) {

    double dx = detX - srcX;
    double dy = detY - srcY;
    double dz = detZ - srcZ;
    double length = sqrt(dx*dx + dy*dy + dz*dz);
    dx /= length; dy /= length; dz /= length;

    double integral = 0.0;
    int numSamples = 500;
    double stepSize = length / numSamples;

    for (int i = 0; i < numSamples; i++) {
        double t = i * stepSize;
        double x = srcX + t * dx;
        double y = srcY + t * dy;
        double z = srcZ + t * dz;
        integral += phantom.getDensity(x, y, z) * stepSize;
    }
    return integral;
}

void FDKReconstructor::cosineWeight(std::vector<std::vector<double>>& projection) {
    double minWeight = 1e10, maxWeight = -1e10;

    for (int row = 0; row < params_.numDetectorRows; row++) {
        for (int col = 0; col < params_.numDetectorCols; col++) {
            double Y = (col - params_.numDetectorCols/2.0) * params_.detectorSpacing;
            double Z = (row - params_.numDetectorRows/2.0) * params_.detectorSpacing;
            double weight = params_.D / sqrt(params_.D * params_.D + Y * Y + Z * Z);

            minWeight = std::min(minWeight, weight);
            maxWeight = std::max(maxWeight, weight);

            projection[row][col] *= weight;
        }
    }

    static bool firstTime = true;
    if (firstTime) {
        std::cout << "\nCosine weighting range: [" << minWeight << ", " << maxWeight << "]" << std::endl;
        firstTime = false;
    }
}

void fft1d(std::vector<std::complex<double>>& data, const bool inverse) {
    int N = data.size();
    auto* in = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * N));
    auto* out = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * N));

    for (int i = 0; i < N; i++) {
        in[i][0] = data[i].real();
        in[i][1] = data[i].imag();
    }

    int direction = inverse ? FFTW_BACKWARD : FFTW_FORWARD;
    fftw_plan plan = fftw_plan_dft_1d(N, in, out, direction, FFTW_ESTIMATE);
    fftw_execute(plan);

    double norm = inverse ? (1.0 / N) : 1.0;
    for (int i = 0; i < N; i++) {
        data[i] = std::complex<double>(out[i][0] * norm, out[i][1] * norm);
    }

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);
}

void FDKReconstructor::filterProjection(std::vector<std::vector<double>>& projection) {
    static bool firstTime = true;
    double minBefore = 1e10, maxBefore = -1e10;
    double minAfter = 1e10, maxAfter = -1e10;

    for (int row = 0; row < params_.numDetectorRows; row++) {
        int N = params_.numDetectorCols;
        int paddedN = 1;
        while (paddedN < N) paddedN *= 2;
        paddedN *= 2;

        std::vector<std::complex<double>> rowData(paddedN, 0.0);
        for (int col = 0; col < N; col++) {
            rowData[col] = projection[row][col];
            if (firstTime && row == params_.numDetectorRows/2) {
                minBefore = std::min(minBefore, projection[row][col]);
                maxBefore = std::max(maxBefore, projection[row][col]);
            }
        }

        fft1d(rowData, false);

        double deltaS = params_.detectorSpacing;
        for (int i = 0; i < paddedN; i++) {
            int k = (i <= paddedN/2) ? i : (i - paddedN);
            double freq = k / (paddedN * deltaS);
            double H = std::abs(freq);
            double fNyq = 1.0 / (2.0 * deltaS);
            if (std::abs(freq) > fNyq) H = 0.0;
            rowData[i] *= H;
        }

        fft1d(rowData, true);

        for (int col = 0; col < N; col++) {
            projection[row][col] = rowData[col].real();
            if (firstTime && row == params_.numDetectorRows/2) {
                minAfter = std::min(minAfter, projection[row][col]);
                maxAfter = std::max(maxAfter, projection[row][col]);
            }
        }
    }

    if (firstTime) {
        std::cout << "\nFilter effect on center row:" << std::endl;
        std::cout << "  Before: [" << minBefore << ", " << maxBefore << "]" << std::endl;
        std::cout << "  After:  [" << minAfter << ", " << maxAfter << "]" << std::endl;
        firstTime = false;
    }
}

void FDKReconstructor::backproject(
    const std::vector<std::vector<std::vector<double>>>& filteredProjections,
    std::vector<std::vector<std::vector<double>>>& volume) const {

    std::cout << "\n=== Starting Backprojection ===" << std::endl;

    int nx = volume[0][0].size();
    int ny = volume[0].size();
    int nz = volume.size();

    std::cout << "Volume size: " << nz << " x " << ny << " x " << nx << std::endl;

    double fov = 40.0;
    double gridSpacing = fov / (nx - 1);
    std::cout << "FOV: " << fov << ", Grid spacing: " << gridSpacing << std::endl;

    int totalVoxels = 0;
    int voxelsWithData = 0;
    double minVal = 1e10, maxVal = -1e10;

    for (int iz = 0; iz < nz; iz++) {
        double z = (iz - (nz-1)/2.0) * gridSpacing;

        for (int iy = 0; iy < ny; iy++) {
            double y = (iy - (ny-1)/2.0) * gridSpacing;

            for (int ix = 0; ix < nx; ix++) {
                double x = (ix - (nx-1)/2.0) * gridSpacing;
                double sum = 0.0;
                int hitCount = 0;

                for (int angleIdx = 0; angleIdx < params_.numAngles; angleIdx++) {
                    double phi = 2.0 * PI * angleIdx / params_.numAngles;
                    double srcX = params_.d * cos(phi);
                    double srcY = params_.d * sin(phi);
                    double srcZ = 0.0;

                    double rx = x - srcX;
                    double ry = y - srcY;
                    double rz = z - srcZ;
                    double rayLength = sqrt(rx*rx + ry*ry + rz*rz);

                    double rdx = rx / rayLength;
                    double rdy = ry / rayLength;
                    double rdz = rz / rayLength;

                    double detDist = params_.D - params_.d;
                    double detCX = -detDist * cos(phi);
                    double detCY = -detDist * sin(phi);
                    double detCZ = 0.0;

                    double normX = cos(phi);
                    double normY = sin(phi);
                    double normZ = 0.0;

                    double denom = rdx * normX + rdy * normY + rdz * normZ;
                    if (fabs(denom) < 1e-10) continue;

                    double t = ((detCX - srcX) * normX +
                               (detCY - srcY) * normY +
                               (detCZ - srcZ) * normZ) / denom;

                    if (t <= 0) continue;

                    double intX = srcX + t * rdx;
                    double intY = srcY + t * rdy;
                    double intZ = srcZ + t * rdz;

                    double dx_det = intX - detCX;
                    double dy_det = intY - detCY;
                    double dz_det = intZ - detCZ;

                    double uX = -sin(phi);
                    double uY = cos(phi);

                    double u = dx_det * uX + dy_det * uY;
                    double v = dz_det;

                    double colFloat = u / params_.detectorSpacing + params_.numDetectorCols / 2.0;
                    double rowFloat = v / params_.detectorSpacing + params_.numDetectorRows / 2.0;

                    int col0 = (int)floor(colFloat);
                    int row0 = (int)floor(rowFloat);

                    if (col0 >= 0 && col0 < params_.numDetectorCols - 1 &&
                        row0 >= 0 && row0 < params_.numDetectorRows - 1) {

                        double dc = colFloat - col0;
                        double dr = rowFloat - row0;

                        double val = (1-dr) * (1-dc) * filteredProjections[angleIdx][row0][col0]
                                   + (1-dr) * dc * filteredProjections[angleIdx][row0][col0+1]
                                   + dr * (1-dc) * filteredProjections[angleIdx][row0+1][col0]
                                   + dr * dc * filteredProjections[angleIdx][row0+1][col0+1];

                        // After computing t and checking t > 0:
                        double L = t; // distance from source to detector intersection
                        double weight = (params_.d * params_.d) / (L * L);
                        sum += weight * val;


                        sum += weight * val;
                        hitCount++;
                    }
                }

                volume[iz][iy][ix] = sum / params_.numAngles;
                totalVoxels++;
                if (hitCount > 0) voxelsWithData++;
                minVal = std::min(minVal, volume[iz][iy][ix]);
                maxVal = std::max(maxVal, volume[iz][iy][ix]);
            }
        }
    }

    std::cout << "Voxels with data: " << voxelsWithData << " / " << totalVoxels << std::endl;
    std::cout << "Reconstruction value range: [" << minVal << ", " << maxVal << "]" << std::endl;
}

std::vector<std::vector<std::vector<double>>>
FDKReconstructor::reconstruct(
    const std::vector<std::vector<std::vector<double>>>& projections) {

    std::cout << "\n=== Starting Reconstruction ===" << std::endl;
    auto filteredProj = projections;

    for (int angleIdx = 0; angleIdx < params_.numAngles; angleIdx++) {
        cosineWeight(filteredProj[angleIdx]);
        filterProjection(filteredProj[angleIdx]);
    }

    int volSize = 99;
    std::vector<std::vector<std::vector<double>>> volume(
        49,
        std::vector<std::vector<double>>(volSize, std::vector<double>(volSize, 0.0))
    );

    backproject(filteredProj, volume);
    return volume;
}