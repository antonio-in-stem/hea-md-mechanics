# LAMMPS compression protocol

## Requirements

Use a LAMMPS executable with the MEAM package and the Huang Hf–Nb–Ta–Ti–Zr parameter files:

```text
library.meam
HfNbTaTiZr.meam
```

The library selection and atom-type mapping both use `Hf Nb Ta Ti Zr`. The protocol uses second-neighbor interactions (`nn2=1` for all 15 element pairs). [Model and command references](../docs/README.md#4-model-references-and-sources).

## Prepare and run

From the repository root, supply the parameter directory and choose a new run directory:

```bash
python -m hea_md stage eq --potential-dir /path/to/parameters --destination runs/eq
cd runs/eq
lmp -in in.compression.lmp -log log.lammps
```

On Windows, substitute the local executable name and parameter path. `stage` copies the selected input and both parameter files, creates `outputs/`, and records their hashes in `run-inputs.json`. It checks library element names and second-neighbor flags before creating the directory. The command prepares files; the explicit `lmp` command starts the calculation.

The trajectory dumps can occupy several gigabytes. Keep each composition in its own run directory. Capture the executable version, per-species counts, log and parameter hashes with the run output.

## Cell, preparation and loading

The cell contains 6750 initially BCC sites in a `51 × 51 × 51 Å` periodic box, aligned with [100], [010] and [001]. `units metal` supplies Å, ps, eV and bar. The timestep is 0.001 ps.

| Preparation stage | Ensemble | Steps | Time (ps) |
|---|---|---:|---:|
| 300 K | NPT | 30000 | 30 |
| 300 to 1500 K | NPT | 300000 | 300 |
| 1500 K | NVT | 30000 | 30 |
| 1500 to 300 K | NPT | 300000 | 300 |
| 300 K | NPT | 300000 | 300 |
| 300 K | NVT | 50000 | 50 |

NPT preparation stages use an isotropic 0-bar target. The timestep-scaled temperature and pressure damping times are 0.1 and 1 ps. The preparation lasts 1.01 ns.

During compression, NPT controls temperature and the lateral x/y pressures. `fix deform` imposes the z rate. Reference lengths are frozen at the beginning of loading. The equiatomic and Zr-28% inputs run 100 ps; the remaining inputs run 200 ps.

## Nominal species assignment

The chemical seed is 1234567. Successive `type/fraction` selections repeat that seed before atoms move, producing nested coordinate-based selections in the referenced implementation. Differences between thresholds define the nominal fractions. Atom counts remain probabilistic, so record the realized populations for a run.

Gaussian velocity seed 12345 is shared across cases. Each composition has one loading record in the numerical dataset. The inputs specify one cell size, loading direction and strain rate.
