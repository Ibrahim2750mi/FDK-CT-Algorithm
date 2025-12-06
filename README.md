## How to Run


### 1. Compile
```bash
cmake -B build -G "MinGW Makefiles" && cmake --build build --config Release
```

### 2. Set parameters

In main.cpp line 200, edit these parameters to fit the test.

```cpp
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
```

### 3a. Run the exe in (CMD)

```bash
build\fdk_recon.exe "projections\proj_%%03d.png" 90
```

### 3b. Run the exe in Powershell

```bash
.\build\fdk_recon.exe "projections\proj_%03d.png" 90
```

The exe expects the projections folder to contain to slices as proj_000.png to proj_089.png. <br>
If eg: the slices are stored in "slices/" and total 360 of them are present then argument will be <br>
`slices/slice_%03d.png 360`. `%03d` means three integer characters. If they are labelled as `slice_89.png` then <br>
change it to `%02d`.