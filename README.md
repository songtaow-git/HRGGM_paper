# HR-GGM

HR-GGM estimates a Gaussian graphical model that augments measured local
features with a small number of prior-guided hidden features and broad global
features. The current implementation is a C++17 project built with CMake and
OpenMP.

This repository is intended to be used as follows:

1. Use the provided MATLAB writer to convert the local-feature-by-sample
   `double` matrix to `Data_whole.bin`, and prepare the remaining files in
   `Data/`.
2. Edit the adjustable parameters at the beginning of `Main.cpp`.
3. Build the project on Linux.
4. Run the `C_GRN` executable and inspect the files in `Result/`.

## Quick start

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 8

cp -r example/Data build/Data
./build/C_GRN 2>&1 | tee build/HRGGM_run.log
```

The final estimated precision matrix is written to:

```text
build/Result/Theta_final.txt
```

## Documentation

See the complete [HR-GGM User Manual](docs/HRGGM_User_Manual.md) for:

- Linux requirements and build instructions;
- exact input requirements and file formats;
- exact formats of all files in `Data/`;
- detailed explanation of C++ zero-based indices;
- parameters that users should edit in `Main.cpp`;
- parameter selection and model-fitting workflow;
- output files and precision-matrix block interpretation;
- troubleshooting guidance.

The MATLAB writer performs only the conversion to the binary format required by
C++. Dataset-specific preprocessing and simulation-data generation should be
provided separately with the corresponding analyses.

## Citation

Please add the final HR-GGM paper citation here before public release.

## License

Please add the selected software license and retain the license notice for the
bundled Eigen dependency before public release.
