## How to Run


### 1. Compile
```bash
cmake -B build -G "MinGW Makefiles" && cmake --build build --config Release
```

### 2. Set parameters

In main.cpp line 58, edit these parameters to fit the test.

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

### 3. Add Path

Edit main.cpp line 74, to have the base path for example if the slices are in `folder/slice_012.png` enter
`auto projections = loadProjectionsOpenCV(params, "folder/slice_");`
IF theres no zero(i.e `slice_12.png`) then change the `setfill` in line 35 to 2.

### 4a. Run the exe in (CMD)

```bash
build\fdk_recon.exe
```

### 4b. Run the exe in Powershell

```bash
.\build\fdk_recon.exe
```