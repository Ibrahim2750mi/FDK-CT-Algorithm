#include "../include/fdk.h"
#include <cmath>
#include <complex>
#include <iostream>

static const double PI = 3.141592653589793;

// ============================================================================
// Constructor
// ============================================================================
FDKReconstructor::FDKReconstructor(const GeometryParams& params)
    : params_(params)
{
    std::cout << "FDK Reconstructor initialized:\n";
    std::cout << "  SID (d): " << params_.d << " mm\n";
    std::cout << "  SDD: " << params_.sdd << " mm\n";
    std::cout << "  Pitch per rotation = " << params_.helicalPitch << " mm\n";
    std::cout << "  Detector = " << params_.numDetectorCols
              << " x " << params_.numDetectorRows << "\n";
}

// ============================================================================
// Compute line integral through phantom (ray caster)
// ============================================================================
double FDKReconstructor::computeLineIntegral(const Phantom& phantom,
                                             double sx, double sy, double sz,
                                             double dx, double dy, double dz) const
{
    double vx = dx - sx;
    double vy = dy - sy;
    double vz = dz - sz;

    double L = std::sqrt(vx*vx + vy*vy + vz*vz);
    if (L < 1e-9) return 0;

    vx /= L; vy /= L; vz /= L;

    double maxDist = L;
    int steps = 800;
    double dt = maxDist / steps;

    double sum = 0.0;

    for (int i = 0; i < steps; i++)
    {
        double t = i * dt;

        double x = sx + t * vx;
        double y = sy + t * vy;
        double z = sz + t * vz;

        // Convert physical mm into normalized [-1,1] region
        double xn = x / 100.0;
        double yn = y / 100.0;
        double zn = z / 100.0;

        double density = phantom.getDensity(xn, yn, zn);

        sum += density * dt;
    }

    return sum;
}

// ============================================================================
// Generate helical projections - FULLY FIXED
// ============================================================================
std::vector<std::vector<std::vector<double>>>
FDKReconstructor::generateProjections(const Phantom& phantom)
{
    std::cout << "Generating helical projections...\n";

    int Nu = params_.numDetectorCols;
    int Nv = params_.numDetectorRows;
    int NA = params_.numAngles;

    double du = params_.detectorSpacing;
    double dv = params_.detectorSpacing;

    double u0 = (Nu - 1) * 0.5;
    double v0 = (Nv - 1) * 0.5;

    std::vector<std::vector<std::vector<double>>>
        proj(NA, std::vector<std::vector<double>>(Nv, std::vector<double>(Nu)));

    for (int ia = 0; ia < NA; ia++)
    {
        double lam = 2 * PI * ia / NA;

        // Source on helical trajectory
        double sx = params_.d * cos(lam);
        double sy = params_.d * sin(lam);
        double sz = params_.h * lam;

        // Direction from source to isocenter (opposite of source position for circular part)
        // For helical cone-beam, we point toward the isocenter plane at this z-height
        double iso_x = 0.0;
        double iso_y = 0.0;
        double iso_z = sz;  // Isocenter at same z as source for this angle

        // Vector from source to isocenter
        double cx_dir = iso_x - sx;
        double cy_dir = iso_y - sy;
        double cz_dir = iso_z - sz;
        double c_norm = sqrt(cx_dir*cx_dir + cy_dir*cy_dir + cz_dir*cz_dir);
        cx_dir /= c_norm;
        cy_dir /= c_norm;
        cz_dir /= c_norm;

        // Detector center at distance sdd from source in direction toward isocenter
        double cx = sx + params_.sdd * cx_dir;
        double cy = sy + params_.sdd * cy_dir;
        double cz = sz + params_.sdd * cz_dir;

        // Detector u-axis (horizontal, perpendicular to radial direction)
        // Points in the tangent direction of the circular trajectory
        double ux = -sin(lam);
        double uy =  cos(lam);
        double uz =  0.0;

        // Detector v-axis (vertical, pointing up in z)
        double vx = 0.0;
        double vy = 0.0;
        double vz = 1.0;

        for (int r = 0; r < Nv; r++)
        {
            double v_mm = (r - v0) * dv;

            for (int c = 0; c < Nu; c++)
            {
                double u_mm = (c - u0) * du;

                // Detector pixel position
                double dx = cx + u_mm * ux + v_mm * vx;
                double dy = cy + u_mm * uy + v_mm * vy;
                double dz = cz + u_mm * uz + v_mm * vz;

                proj[ia][r][c] = computeLineIntegral(phantom, sx, sy, sz, dx, dy, dz);
            }
        }

        if ((ia+1) % 10 == 0)
            std::cout << "  " << (ia+1) << "/" << NA << " projections\n";
    }

    return proj;
}

// ============================================================================
// Cosine weighting (required by FDK)
// ============================================================================
void FDKReconstructor::cosineWeight(std::vector<std::vector<double>>& P)
{
    int Nu = params_.numDetectorCols;
    int Nv = params_.numDetectorRows;
    double du = params_.detectorSpacing;

    double u0 = (Nu - 1) * 0.5;
    double v0 = (Nv - 1) * 0.5;

    for (int r = 0; r < Nv; r++)
    {
        double v = (r - v0) * du;

        for (int c = 0; c < Nu; c++)
        {
            double u = (c - u0) * du;

            double w = params_.d / sqrt(params_.d * params_.d + u*u + v*v);
            P[r][c] *= w;
        }
    }
}

// ============================================================================
// Ramp filter via FFT
// ============================================================================
void fft1d(std::vector<std::complex<double>>& data, const bool inverse) {
    const std::size_t N = data.size();
    if (N == 0) return;

    // --- (Optional) check: N should be a power of 2 for this implementation ---
    if ((N & (N - 1)) != 0) {
        // Not a power of two: you can either handle it with a slow O(N^2) DFT
        // or assert/throw. For now we'll just do a naive DFT fallback.
        std::vector<std::complex<double>> out(N);
        const double sign = inverse ? 1.0 : -1.0;
        const double PI = 3.14159265358979323846;

        for (std::size_t k = 0; k < N; ++k) {
            std::complex<double> sum(0.0, 0.0);
            for (std::size_t n = 0; n < N; ++n) {
                double angle = 2.0 * PI * k * n / static_cast<double>(N);
                std::complex<double> w(std::cos(sign * angle),
                                       std::sin(sign * angle));
                sum += data[n] * w;
            }
            out[k] = sum;
        }

        double norm = inverse ? (1.0 / static_cast<double>(N)) : 1.0;
        for (std::size_t k = 0; k < N; ++k) {
            data[k] = out[k] * norm;
        }
        return;
    }

    // --- Bit-reversal permutation ---
    std::size_t logN = 0;
    while ((1u << logN) < N) ++logN;

    for (std::size_t i = 1, j = 0; i < N; ++i) {
        std::size_t bit = N >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) {
            std::swap(data[i], data[j]);
        }
    }

    const double PI = 3.14159265358979323846;
    // FFTW convention: forward uses exp(-2πikn/N), backward uses exp(+2πikn/N)
    const double sign = inverse ? 1.0 : -1.0;

    // --- Iterative Cooley–Tukey FFT ---
    for (std::size_t len = 2; len <= N; len <<= 1) {
        double angle = 2.0 * PI / static_cast<double>(len) * sign;
        std::complex<double> wlen(std::cos(angle), std::sin(angle));

        for (std::size_t i = 0; i < N; i += len) {
            std::complex<double> w(1.0, 0.0);
            std::size_t half = len >> 1;
            for (std::size_t j = 0; j < half; ++j) {
                std::complex<double> u = data[i + j];
                std::complex<double> v = data[i + j + half] * w;
                data[i + j]         = u + v;
                data[i + j + half]  = u - v;
                w *= wlen;
            }
        }
    }

    // --- Normalization (match your FFTW wrapper) ---
    if (inverse) {
        double norm = 1.0 / static_cast<double>(N);
        for (auto& x : data) {
            x *= norm;
        }
    }
}



std::vector<double> FDKReconstructor::rampFilter(const std::vector<double>& sig) const
{
    int N = static_cast<int>(sig.size());
    if (N == 0) return {};

    // M is next power of two >= 2N
    int M = 1;
    while (M < 2 * N) M <<= 1;

    // Prepare complex buffer with zero-padding
    std::vector<std::complex<double>> buf(M);
    for (int i = 0; i < N; ++i) {
        buf[i] = std::complex<double>(sig[i], 0.0);
    }
    for (int i = N; i < M; ++i) {
        buf[i] = std::complex<double>(0.0, 0.0);
    }

    // Forward FFT (no scaling inside fft1d for forward)
    fft1d(buf, /*inverse=*/false);

    // Apply ramp filter in frequency domain
    for (int k = 0; k < M; ++k) {
        double omega;
        if (k <= M / 2) {
            omega = 2.0 * PI * k / (M * params_.detectorSpacing);
        } else {
            omega = 2.0 * PI * (k - M) / (M * params_.detectorSpacing);
        }

        double H = std::abs(omega);  // |ω| ramp
        buf[k] *= H;
    }

    // Inverse FFT: our fft1d *already* divides by M
    fft1d(buf, /*inverse=*/true);

    // Extract real part (no extra /M, because inverse already normalized)
    std::vector<double> filtered(N);
    for (int i = 0; i < N; ++i) {
        filtered[i] = buf[i].real();
    }

    return filtered;
}

// ============================================================================
// Filtering (row-wise ramp filter)
// ============================================================================
void FDKReconstructor::filterProjection(std::vector<std::vector<double>>& P)
{
    int Nv = params_.numDetectorRows;
    int Nu = params_.numDetectorCols;

    std::vector<std::vector<double>> out(Nv, std::vector<double>(Nu));

    for (int r = 0; r < Nv; r++)
    {
        std::vector<double> row(Nu);
        for (int c = 0; c < Nu; c++)
            row[c] = P[r][c];

        std::vector<double> fr = rampFilter(row);

        for (int c = 0; c < Nu; c++)
            out[r][c] = fr[c];
    }

    P.swap(out);
}

// ============================================================================
// Backprojection - FULLY FIXED
// ============================================================================
void FDKReconstructor::backproject(
    const std::vector<std::vector<std::vector<double>>>& FP,
    std::vector<std::vector<std::vector<double>>>& vol) const
{
    int N = params_.reconSize;
    double s = params_.reconSpacing;

    int Nu = params_.numDetectorCols;
    int Nv = params_.numDetectorRows;

    double u0 = (Nu - 1) * 0.5;
    double v0 = (Nv - 1) * 0.5;
    double du = params_.detectorSpacing;

    int NA = params_.numAngles;

    for (int iz = 0; iz < N; iz++)
    {
        double z = (iz - N/2.0) * s;

        for (int iy = 0; iy < N; iy++)
        {
            double y = (iy - N/2.0) * s;

            for (int ix = 0; ix < N; ix++)
            {
                double x = (ix - N/2.0) * s;

                double sum = 0;

                for (int ia = 0; ia < NA; ia++)
                {
                    double lam = 2*PI * ia / NA;

                    double sx = params_.d * cos(lam);
                    double sy = params_.d * sin(lam);
                    double sz = params_.h * lam;

                    // Vector from source to voxel
                    double dx = x - sx;
                    double dy = y - sy;
                    double dz = z - sz;

                    // Direction from source to isocenter
                    double iso_x = 0.0;
                    double iso_y = 0.0;
                    double iso_z = sz;

                    double cx_dir = iso_x - sx;
                    double cy_dir = iso_y - sy;
                    double cz_dir = iso_z - sz;
                    double c_norm = sqrt(cx_dir*cx_dir + cy_dir*cy_dir + cz_dir*cz_dir);
                    cx_dir /= c_norm;
                    cy_dir /= c_norm;
                    cz_dir /= c_norm;

                    // Distance along central ray direction
                    double denom = dx * cx_dir + dy * cy_dir + dz * cz_dir;

                    if (denom < 1e-6) continue;  // Behind source or too close

                    // Detector axes (same as forward projection)
                    double ux = -sin(lam);
                    double uy =  cos(lam);
                    double uz =  0.0;

                    double vx = 0.0;
                    double vy = 0.0;
                    double vz = 1.0;

                    // Project onto detector coordinate system
                    double U = (params_.sdd / denom) * (dx * ux + dy * uy + dz * uz);
                    double V = (params_.sdd / denom) * (dx * vx + dy * vy + dz * vz);

                    double cu = U / du + u0;
                    double cv = V / du + v0;

                    if (cu < 0 || cu >= Nu-1 ||
                        cv < 0 || cv >= Nv-1)
                        continue;

                    int c0 = floor(cu);
                    int r0 = floor(cv);
                    double dc = cu - c0;
                    double dr = cv - r0;

                    const auto& P = FP[ia];

                    double p00 = P[r0][c0];
                    double p01 = P[r0][c0+1];
                    double p10 = P[r0+1][c0];
                    double p11 = P[r0+1][c0+1];

                    double val =
                        (1-dc)*(1-dr)*p00 +
                         dc   *(1-dr)*p01 +
                        (1-dc)*   dr*p10 +
                         dc *     dr*p11;

                    double weight = (params_.d * params_.d)/(denom*denom);

                    sum += weight * val;
                }

                vol[iz][iy][ix] = sum * (2*PI / NA);
            }
        }

        std::cout << "Slice " << iz+1 << "/" << N << "\n";
    }
}

// ============================================================================
// Full reconstruction
// ============================================================================
std::vector<std::vector<std::vector<double>>>
FDKReconstructor::reconstruct(const std::vector<std::vector<std::vector<double>>>& projections)
{
    auto P = projections;

    std::cout << "Applying cosine weighting...\n";
    for (auto& pr : P)
        cosineWeight(pr);

    std::cout << "Filtering...\n";
    for (auto& pr : P)
        filterProjection(pr);

    std::vector<std::vector<std::vector<double>>>
        vol(params_.reconSize,
            std::vector<std::vector<double>>(params_.reconSize,
                std::vector<double>(params_.reconSize, 0.0)));

    std::cout << "Backprojecting...\n";
    backproject(P, vol);

    return vol;
}

// ============================================================================
// Parker window placeholder (not used here)
// ============================================================================
double FDKReconstructor::parkerWindow(double beta, double gamma) const
{
    return 1.0;
}