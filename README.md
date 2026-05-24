# MalmstenIntegralSequences

Wolfram Language package for Malmsten-type integral sequences and nested-sum identities.

This repository contains code accompanying Abdulhafeez A. Abdulsalam's work on Malmsten's integral and related integral sequences. The package collects reusable symbolic computations for closed forms, coefficient vectors, nested sums, and reductions involving Hurwitz-zeta derivative terms.

## Contents

- `MalmstenIntegralSequences.zip` contains the package source.
- The package includes kernel code, tests, and a short documentation guide.

## Main Functions

- `SechPowerIntegral[n]`
- `MalmstenSequence[n, a, b]`
- `ChiSequence[n]`
- `FSequence[n]`
- `LambdaSequence[n]`
- `DeltaSequence[n]`
- `NestedChiCoefficients[N, lows]`
- `NestedFRHS[N, lows]`
- `NestedFSymbolicRHS[N, lows]`
- `EvaluateFRHS[expr]`
- `NestedChiIntegrandPolynomial[N, lows, x]`
- `ReduceQuarterZetaTerms[expr]`

## Basic Usage

After downloading and unzipping the package:

```wl
SetDirectory["/path/to/the/unzipped/folder"];

Get["MalmstenIntegralSequences/Kernel/MalmstenIntegralSequences.wl"];

```

Example:

```wl
MalmstenSequence[1]

NestedFSymbolicRHS[5, {2, 2, 2, 2, 2}]

EvaluateFRHS[NestedFSymbolicRHS[5, {2, 2, 2, 2, 2}]]
```

## Verification

To run the notebook-friendly package check:

```wl
Get["MalmstenIntegralSequences/Tests/VerifyPackage.wls"]
```

A successful run should print a `PASS` message or show a test report with 100% success.

## Citation

If you use this package, please cite:

Abdulhafeez A. Abdulsalam *MalmstenIntegralSequences: A Wolfram Language package for Malmsten-type integral sequences*, version 1.0, 2026.

A formal DOI citation will be added after archival release.

## License

MIT License.
# MalmstenIntegralSequences
Initial package release
