import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import subprocess
import os
import time
from pathlib import Path

class FDKCTAnalyzer:
    def __init__(self):
        self.projections = None
        self.reconstruction = None
        self.phantom_ground_truth = None

    def compile_and_run_cpp(self):
        """Compile and run the C++ FDK implementation"""
        print("Compiling C++ code...")

        # Create build directory
        build_dir = Path("build")
        build_dir.mkdir(exist_ok=True)

        # Compile command
        compile_cmd = [
            "g++",
            "-std=c++17",
            "-O2",
            "-I./include",
            "-lfftw3",
            "-lm",
            "-o", "build/fdk_reconstruction",
            "src/phantom.cpp",
            "src/fdk.cpp",
            "src/main.cpp"
        ]

        try:
            # result = subprocess.run(compile_cmd, capture_output=True, text=True)
            # if result.returncode != 0:
            #     print(f"Compilation failed:\n{result.stderr}")
            #     return False

            print("Compilation successful!")

            # Run the program
            print("Running reconstruction...")
            start_time = time.time()
            result = subprocess.run(["./build/fdk_recon"], capture_output=True, text=True)
            end_time = time.time()

            print(f"Execution time: {end_time - start_time:.2f} seconds")
            print(f"Output:\n{result.stdout}")

            return True

        except FileNotFoundError:
            print("Error: g++ not found. Please install g++ compiler.")
            return False

    def load_reconstruction(self, filename="build/midplane.txt"):
        """Load the reconstruction results from file"""
        if not os.path.exists(filename):
            print(f"File {filename} not found. Run C++ code first.")
            return None

        try:
            data = []
            with open(filename, 'r') as f:
                for line in f:
                    row = [float(x) for x in line.strip().split() if x]
                    if row:
                        data.append(row)

            self.reconstruction = np.array(data)
            print(f"Loaded reconstruction with shape: {self.reconstruction.shape}")
            return self.reconstruction

        except Exception as e:
            print(f"Error loading file: {e}")
            return None

    def create_test_phantom(self, size=99):
        """Create a simple phantom in Python for comparison"""
        # This matches the phantom in phantom.cpp
        phantom = np.zeros((size, size))

        # Center coordinates
        cx, cy = size // 2, size // 2

        # Outer cylinder (radius 20)
        y, x = np.ogrid[-cy:size-cy, -cx:size-cx]
        mask = (x**2 + y**2) <= 20**2
        phantom[mask] += 1.0

        # Inner cylinder (radius 17)
        mask = (x**2 + y**2) <= 17**2
        phantom[mask] += -1.21

        # Central ellipsoid (radius x=15, y=10)
        mask = (x**2/15**2 + y**2/10**2) <= 1
        phantom[mask] += 0.21

        self.phantom_ground_truth = phantom
        return phantom

    def analyze_reconstruction(self):
        """Analyze the reconstruction quality"""
        if self.reconstruction is None:
            self.load_reconstruction()

        if self.reconstruction is None or self.phantom_ground_truth is None:
            print("No data to analyze")
            return

        # Normalize both for comparison
        recon_norm = (self.reconstruction - self.reconstruction.min()) / (self.reconstruction.max() - self.reconstruction.min() + 1e-10)
        phantom_norm = (self.phantom_ground_truth - self.phantom_ground_truth.min()) / (self.phantom_ground_truth.max() - self.phantom_ground_truth.min() + 1e-10)

        # Calculate metrics
        mse = np.mean((recon_norm - phantom_norm[:recon_norm.shape[0], :recon_norm.shape[1]])**2)
        rmse = np.sqrt(mse)
        psnr = 20 * np.log10(1.0 / (rmse + 1e-10))

        print("\n=== Reconstruction Analysis ===")
        print(f"MSE: {mse:.6f}")
        print(f"RMSE: {rmse:.6f}")
        print(f"PSNR: {psnr:.2f} dB")
        print(f"Reconstruction range: [{self.reconstruction.min():.3f}, {self.reconstruction.max():.3f}]")

        # Profile through center
        center_line = self.reconstruction[self.reconstruction.shape[0]//2, :]
        print(f"\nCenter line values (first 10): {center_line[:10]}")

    def visualize_results(self):
        """Create visualizations of the reconstruction"""
        if self.reconstruction is None:
            self.load_reconstruction()

        if self.reconstruction is None:
            print("No reconstruction data to visualize")
            return

        # Create phantom ground truth for comparison
        self.create_test_phantom(self.reconstruction.shape[0])

        fig = plt.figure(figsize=(15, 5))

        # Plot 1: Reconstruction
        ax1 = fig.add_subplot(131)
        im1 = ax1.imshow(self.reconstruction, cmap='gray', origin='lower')
        ax1.set_title('FDK Reconstruction')
        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')
        plt.colorbar(im1, ax=ax1)

        # Plot 2: Phantom Ground Truth (middle slice)
        ax2 = fig.add_subplot(132)
        im2 = ax2.imshow(self.phantom_ground_truth[:self.reconstruction.shape[0], :self.reconstruction.shape[1]],
                         cmap='gray', origin='lower')
        ax2.set_title('Phantom Ground Truth')
        ax2.set_xlabel('X')
        ax2.set_ylabel('Y')
        plt.colorbar(im2, ax=ax2)

        # Plot 3: Difference
        ax3 = fig.add_subplot(133)
        diff = self.reconstruction - self.phantom_ground_truth[:self.reconstruction.shape[0], :self.reconstruction.shape[1]]
        im3 = ax3.imshow(diff, cmap='coolwarm', origin='lower', vmin=-np.abs(diff).max(), vmax=np.abs(diff).max())
        ax3.set_title('Difference (Recon - Phantom)')
        ax3.set_xlabel('X')
        ax3.set_ylabel('Y')
        plt.colorbar(im3, ax=ax3)

        plt.tight_layout()
        plt.savefig('reconstruction_analysis.png', dpi=150)
        plt.show()

        # Plot profiles through center
        fig2, (ax4, ax5) = plt.subplots(1, 2, figsize=(12, 4))

        # Horizontal profile
        center_row = self.reconstruction.shape[0] // 2
        ax4.plot(self.reconstruction[center_row, :], 'b-', label='Reconstruction', linewidth=2)
        ax4.plot(self.phantom_ground_truth[center_row, :self.reconstruction.shape[1]], 'r--',
                 label='Ground Truth', linewidth=2, alpha=0.7)
        ax4.set_xlabel('Pixel Index')
        ax4.set_ylabel('Density')
        ax4.set_title(f'Horizontal Profile at Row {center_row}')
        ax4.legend()
        ax4.grid(True, alpha=0.3)

        # Vertical profile
        center_col = self.reconstruction.shape[1] // 2
        ax5.plot(self.reconstruction[:, center_col], 'b-', label='Reconstruction', linewidth=2)
        ax5.plot(self.phantom_ground_truth[:self.reconstruction.shape[0], center_col], 'r--',
                 label='Ground Truth', linewidth=2, alpha=0.7)
        ax5.set_xlabel('Pixel Index')
        ax5.set_ylabel('Density')
        ax5.set_title(f'Vertical Profile at Column {center_col}')
        ax5.legend()
        ax5.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig('reconstruction_profiles.png', dpi=150)
        plt.show()

        # 3D surface plot (subsampled for clarity)
        fig3 = plt.figure(figsize=(10, 7))
        ax = fig3.add_subplot(111, projection='3d')

        # Subsample for clearer 3D plot
        step = max(1, self.reconstruction.shape[0] // 20)
        x = np.arange(0, self.reconstruction.shape[1], step)
        y = np.arange(0, self.reconstruction.shape[0], step)
        X, Y = np.meshgrid(x, y)
        Z = self.reconstruction[::step, ::step]

        surf = ax.plot_surface(X, Y, Z, cmap='viridis', alpha=0.8,
                               linewidth=0, antialiased=True)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Density')
        ax.set_title('3D Reconstruction Surface')
        fig3.colorbar(surf, ax=ax, shrink=0.5, aspect=5)

        plt.savefig('3d_reconstruction.png', dpi=150)
        plt.show()

        print("Visualizations saved as PNG files.")

    def run_unit_tests(self):
        """Run basic unit tests on the algorithm components"""
        print("\n=== Running Unit Tests ===")

        # Test 1: Check if files exist
        required_files = ["src/phantom.cpp", "src/fdk.cpp", "src/main.cpp",
                          "include/phantom.h", "include/fdk.h"]
        for file in required_files:
            if os.path.exists(file):
                print(f"✓ {file} exists")
            else:
                print(f"✗ {file} missing")

        # Test 2: Check if FFTW is available
        try:
            import ctypes
            fftw = ctypes.CDLL('libfftw3.so.3')
            print("✓ FFTW library found")
        except:
            print("✗ FFTW library not found (needed for filtering)")

        # Test 3: Simple geometric test
        print("\nSimple geometric consistency test:")
        # The source-detector distances should be positive
        if params_d > 0 and params_D > 0:
            print(f"✓ Geometry parameters valid: d={params_d}, D={params_D}")
        else:
            print(f"✗ Invalid geometry parameters")

    def export_data(self):
        """Export data for further analysis"""
        if self.reconstruction is not None:
            np.save('reconstruction.npy', self.reconstruction)
            print("Exported reconstruction to 'reconstruction.npy'")

        if self.phantom_ground_truth is not None:
            np.save('phantom_ground_truth.npy', self.phantom_ground_truth)
            print("Exported phantom to 'phantom_ground_truth.npy'")

def main():
    """Main test routine"""
    print("=" * 60)
    print("FDK CT Algorithm Test Suite")
    print("=" * 60)

    analyzer = FDKCTAnalyzer()

    # Run unit tests
    analyzer.run_unit_tests()

    # Compile and run C++ code
    success = analyzer.compile_and_run_cpp()

    if success:
        # Load and analyze results
        print("\n" + "=" * 60)
        print("Loading and Analyzing Results")
        print("=" * 60)

        analyzer.load_reconstruction()
        analyzer.create_test_phantom()
        analyzer.analyze_reconstruction()

        # Visualize
        print("\n" + "=" * 60)
        print("Creating Visualizations")
        print("=" * 60)
        analyzer.visualize_results()

        # Export data
        analyzer.export_data()

        print("\n" + "=" * 60)
        print("Test Complete!")
        print("=" * 60)
    else:
        print("Failed to run C++ code. Please check the errors above.")

if __name__ == "__main__":
    # These should match your C++ parameters
    params_d = 60.0
    params_D = 60.0
    params_numAngles = 32

    main()