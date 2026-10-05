# D4 Specification - Group-I/Group-II chloride initial melt-cell builder

**Authority domain:** D4 public software contract  
**Authority state:** accepted-current reconstruction for executable version 1.3.4  
**Governing protocol:** SSDP 6.6.0  
**Parents:** `SCIENTIFIC_METHOD.md`, `NUMERICAL_METHOD.md`, `ARCHITECTURE.md`

## 1. Entry point and version

The supported executable entry point is:

```bash
python build_mixed_alkali_chloride_rocksalt.py ...
```

`--version` reports executable version `1.3.4`.

## 2. Required cell and temperature inputs

`--nx`, `--ny`, and `--nz` are required positive integers. The cation-site count is `4*nx*ny*nz`.

Exactly one of `--temp-k` or `--temp-c` is required. Celsius/Kelvin conversion is `T_K = T_C + 273.15`. Absolute temperature must be positive in Kelvin.

## 3. Composition input

A run chooses exactly one composition mode:

- mole fractions: `--x-li`, `--x-na`, ..., `--x-ba`; or
- explicit cation counts: `--n-li`, `--n-na`, ..., `--n-ba`.

Fractions must be finite and nonnegative. Counts must be nonnegative integers.

### Legacy Li/Na/K completion

For compositions using only Li/Na/K, the historical convenience behavior remains:

- one supplied fraction must equal 1;
- two supplied fractions complete the third to 1;
- three supplied fractions must sum to 1;
- analogous rules apply to explicit counts with total `N`.

### Extended compositions

If any of Rb, Cs, Be, Mg, Ca, Sr, or Ba is supplied, omitted species are zero and all supplied fractions must sum to 1, or all supplied counts must sum to `N`.

Fraction inputs are converted to integer counts by the D2 apportionment rule in `NUMERICAL_METHOD.md`; all reported "achieved fractions" are derived from those counts.

## 4. Thermodynamic options

`--ideal-only` disables only the legacy Li/Na/K binary excess contributions. Pure-component interpolation remains active.

`--hydrostatic-strain` and alias `--volume-strain` accept either:

- fractional volume strain, e.g. `0.05`; or
- percentage form, e.g. `+5%`.

The value must be greater than `-1` (-100%). Default is `0.0`.

## 5. Randomization

`--seed` accepts an integer and defaults to `12345`. It controls cation-site order and, when divalent species are present, selection of auxiliary chloride sites.

## 6. Output selection

`--format` is one of `poscar`, `cif`, `lammps`, or `all`; default is `all`.

`--output-prefix` defaults to `mixed_group12_chloride`. Unless `--no-auto-temp-tag` is set, a temperature suffix is appended (or a previous compatible suffix is replaced).

For LAMMPS output, `--lammps-atom-style` is `atomic` or `charge`, default `atomic`.

- `atomic` writes atom ID, type, and Cartesian coordinates.
- `charge` additionally writes per-atom formal charge `+1`, `+2`, or `-1` as appropriate.

LAMMPS users must select the corresponding `atom_style` before `read_data` [R9].

## 7. Required scientific output invariants

For every successful run:

- cation counts sum exactly to `N`;
- chloride count equals `sum_i z_i N_i`;
- formal net ionic charge is zero;
- thermodynamic calculations use achieved fractions;
- finite-cell volume follows `NUMERICAL_METHOD.md`;
- lattice determinant equals the intended total volume within floating-point/text-serialization roundoff;
- generated coordinate count matches the reported species counts.

## 8. Runtime diagnostics

The CLI reports:

- input/composition mode;
- cell dimensions and cation-site count;
- chloride and total atom counts;
- net formal ionic charge;
- target temperature;
- achieved cation counts/fractions;
- equilibrium and strained molar volumes/densities;
- cell volumes and lattice dimensions;
- excess-model state;
- topology state;
- LAMMPS atom style when applicable;
- empirical validity warnings;
- written paths.

Empirical validity warnings are nonfatal and observational. Under `--ideal-only`, legacy binary excess-model warnings must be absent while applicable pure-component warnings remain.

## 9. Output-format contract

### POSCAR

The writer emits lattice vectors, species order, species counts, and fractional (`Direct`) coordinates in a VASP-compatible POSCAR representation [R7].

### CIF

The writer emits a periodic P1 CIF containing cell lengths/angles and fractional atom sites [R8].

### LAMMPS data

The writer emits an orthorhombic data file with atom count/types, box bounds, masses, and an `Atoms` section compatible with the selected `atomic` or `charge` style [R9].

## 10. Error contract

The CLI terminates with an actionable error for mutually inconsistent composition modes, invalid cell dimensions, invalid fractions/counts, invalid absolute temperature, invalid hydrostatic strain, unsupported internal state, or failed internal charge/count invariants.

No automatic fallback may silently alter the requested scientific state to force successful output.

## 11. Evidence boundary

`tests/test_builder.py` is executable evidence for this specification and its upstream concretization; it is not the source of scientific truth. In particular, compatibility fixtures protect the frozen Li/Na/K field, while external references govern source-backed empirical claims.

The v1.3.4 suite contains both assembled diagnostic directions:

- active extrapolated legacy binary model -> warning must survive through `main()`;
- inactive binary model under `--ideal-only` -> binary warning must be suppressed.

## References

Exact locators are in `REFERENCES.md`. Serialization requirements use [R7]-[R9]; SI conversion uses [R10].
