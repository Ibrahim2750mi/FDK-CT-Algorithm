#include "../include/fdk.h"
#include <cmath>
#include <algorithm>
#include <iostream>

const double PI = 3.14159265358979323846;

FDKReconstructor::FDKReconstructor(const GeometryParams& params)
    : params_(params) {}

std::vector<std::vector<std::vector<double>>>
FDKReconstructor::generateProjections(const Phantom& phantom) {

    // projections[angle][row][col]
    std::vector<std::vector<std::vector<double>>> projections(
        params_.numAngles,
        std::vector<std::vector<double>>(
            params_.numDetectorRows,
            std::vector<double>(params_.numDetectorCols, 0.0)
        )
    );

    for (int angleIdx = 0; angleIdx < params_.numAngles; angleIdx++) {
        double phi = 2.0 * PI * angleIdx / params_.numAngles;

        // Source position (rotates around z-axis)
        double srcX = params_.d * cos(phi);
        double srcY = params_.d * sin(phi);
        double srcZ = params_.h * phi;

        for (int row = 0; row < params_.numDetectorRows; row++) {
            for (int col = 0; col < params_.numDetectorCols; col++) {

                // Detector coordinates (Section 2, Fig 1)
                double Y = (col - params_.numDetectorCols/2.0) * params_.detectorSpacing;
                double Z = (row - params_.numDetectorRows/2.0) * params_.detectorSpacing;

                // Detector position (opposite side of rotation axis)
                // Detector center radius relative to origin (treat params_.D as source-to-detector)
                double detCenterRadius = params_.D + params_.d; // could be negative (detector past origin)
                double detX = detCenterRadius * cos(phi) + Y * sin(phi);
                double detY = detCenterRadius * sin(phi) - Y * cos(phi);
                double detZ = (params_.h * phi + Z);


                // Compute line integral (Eq. 4)
                projections[angleIdx][row][col] =
                    computeLineIntegral(phantom, srcX, srcY, srcZ,
                                       detX, detY, detZ);
            }
        }
    }

    return projections;
}

double FDKReconstructor::computeLineIntegral(
    const Phantom& phantom,
    double srcX, double srcY, double srcZ,
    double detX, double detY, double detZ) {

    // Ray direction
    double dx = detX - srcX;
    double dy = detY - srcY;
    double dz = detZ - srcZ;
    double length = sqrt(dx*dx + dy*dy + dz*dz);
    dx /= length; dy /= length; dz /= length;

    // Numerical integration along ray
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
    // From before Eq. 31: multiply by d / sqrt(d^2 + Y^2 + Z^2)

    for (int row = 0; row < params_.numDetectorRows; row++) {
        for (int col = 0; col < params_.numDetectorCols; col++) {

            double Y = (col - params_.numDetectorCols/2.0) * params_.detectorSpacing;
            double Z = (row - params_.numDetectorRows/2.0) * params_.detectorSpacing;

            double weight = params_.d / sqrt(params_.d * params_.d +
                                            Y * Y + Z * Z);

            projection[row][col] *= weight;
        }
    }
}


#include <complex>
#include <vector>

#include <fftw3.h>

void fft1d(std::vector<std::complex<double>>& data, const bool inverse) {
    int N = data.size();

    auto* in = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * N));
    auto* out = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * N));

    // Copy data to FFTW format
    for (int i = 0; i < N; i++) {
        in[i][0] = data[i].real();
        in[i][1] = data[i].imag();
    }

    // Create plan and execute
    int direction = inverse ? FFTW_BACKWARD : FFTW_FORWARD;
    fftw_plan plan = fftw_plan_dft_1d(N, in, out, direction, FFTW_ESTIMATE);
    fftw_execute(plan);

    // Copy results back
    double norm = inverse ? (1.0 / N) : 1.0;
    for (int i = 0; i < N; i++) {
        data[i] = std::complex<double>(out[i][0] * norm, out[i][1] * norm);
    }

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);
}

void FDKReconstructor::filterProjection(
    std::vector<std::vector<double>>& projection) {

    // Eq. 32: Ramp filter in frequency domain
    // Filter each row independently (Y direction)

    for (int row = 0; row < params_.numDetectorRows; row++) {
        int N = params_.numDetectorCols;

        // Pad to next power of 2 for FFT efficiency
        int paddedN = 1;
        while (paddedN < N) paddedN *= 2;
        paddedN *= 2; // Extra padding to avoid wraparound

        std::vector<std::complex<double>> rowData(paddedN, 0.0);

        // Copy row data
        for (int col = 0; col < N; col++) {
            rowData[col] = projection[row][col];
        }

        // FFT
        fft1d(rowData, false);

        // Apply ramp filter: |ω| (Eq. 16, 32)
        double omega_max = PI / params_.detectorSpacing;
        for (int i = 0; i < paddedN; i++) {
            double omega = (2.0 * PI * i) / (paddedN * params_.detectorSpacing);
            if (i > paddedN/2) {
                omega = (2.0 * PI * (i - paddedN)) / (paddedN * params_.detectorSpacing);
            }

            // Ramp filter with Shepp-Logan window (mentioned in paper)
            double filter = fabs(omega);
            if (omega != 0) {
                filter *= fabs(sin(omega * params_.detectorSpacing / 2.0) /
                              (omega * params_.detectorSpacing / 2.0));
            }

            // Bandlimit
            if (fabs(omega) > omega_max) {
                filter = 0.0;
            }

            rowData[i] *= filter;
        }

        // Inverse FFT
        fft1d(rowData, true);

        // Copy back
        for (int col = 0; col < N; col++) {
            projection[row][col] = rowData[col].real();
        }
    }
}


void FDKReconstructor::backproject(
    const std::vector<std::vector<std::vector<double>>>& filteredProjections,
    std::vector<std::vector<std::vector<double>>>& volume) const {

    // Eq. 28: Main reconstruction formula

    int nx = volume[0][0].size();
    int ny = volume[0].size();
    int nz = volume.size();

    // Reconstruction grid (assume centered at origin)

    for (int iz = 0; iz < nz; iz++) {
        double gridSpacing = 0.404;
        double z = (iz - nz/2.0) * gridSpacing;

        for (int iy = 0; iy < ny; iy++) {
            double y = (iy - ny/2.0) * gridSpacing;

            for (int ix = 0; ix < nx; ix++) {
                double x = (ix - nx/2.0) * gridSpacing;

                double sum = 0.0;

                // Loop over all projection angles
                for (int angleIdx = 0; angleIdx < params_.numAngles; angleIdx++) {
                    double phi = 2.0 * PI * angleIdx / params_.numAngles;

                    // Source position
                    double srcX = params_.d * cos(phi);
                    double srcY = params_.d * sin(phi);
                    double srcZ = params_.h * phi;

                    // Vector from source to reconstruction point
                    double rx = x - srcX;
                    double ry = y - srcY;
                    double rz = z - srcZ;

                    // Unit vector along x-hat direction (from source to axis)
                    double xhatX = -cos(phi);
                    double xhatY = -sin(phi);

                    // Distance weighting factor
                    double r_dot_xhat = rx * xhatX + ry * xhatY;
                    double U = params_.d + r_dot_xhat;

                    // Check if point is visible from this source position
                    // if (U <= 0) continue;  // Behind the source

                    double weight = (params_.d * params_.d) / (U * U);

                    // Project point onto detector
                    double Y = (params_.d / U) * (rx * sin(phi) - ry * cos(phi));
                    double Z = (params_.d / U) * rz;

                    // Convert to detector indices
                    double colFloat = Y / params_.detectorSpacing + params_.numDetectorCols / 2.0;
                    double rowFloat = Z / params_.detectorSpacing + params_.numDetectorRows / 2.0;

                    // Bilinear interpolation
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

                        sum += weight * val;
                    }
                }

                volume[iz][iy][ix] = sum / (2.0 * PI * params_.numAngles);
            }
        }

        // Progress indicator
        if ((iz + 1) % 10 == 0) {
            std::cout << "  Slice " << (iz + 1) << "/" << nz << std::endl;
        }
    }
}


std::vector<std::vector<std::vector<double>>>
FDKReconstructor::reconstruct(
    const std::vector<std::vector<std::vector<double>>>& projections) {

    // Make a copy for filtering
    auto filteredProj = projections;

    // Process each projection
    for (int angleIdx = 0; angleIdx < params_.numAngles; angleIdx++) {
        // Cosine weighting
        cosineWeight(filteredProj[angleIdx]);

        // Filtering
        filterProjection(filteredProj[angleIdx]);
    }

    // Initialize reconstruction volume
    int volSize = 99; // From paper: 99x99x49
    std::vector<std::vector<std::vector<double>>> volume(
        49,  // z
        std::vector<std::vector<double>>(
            volSize,  // y
            std::vector<double>(volSize, 0.0)  // x
        )
    );

    // Backprojection
    backproject(filteredProj, volume);

    return volume;
}