# Group-I / Group-II chloride initial melt-cell builder

**Executable version:** 1.3.4  
**Governing SSDP:** 6.6.0

This package builds periodic, charge-neutral starting cells for molten chloride simulations using a temperature/composition-dependent volume model. It supports:

- Group I: `Li+`, `Na+`, `K+`, `Rb+`, `Cs+`
- Group II: `Be2+`, `Mg2+`, `Ca2+`, `Sr2+`, `Ba2+`

The generated coordinates are compact starting templates for later melting/equilibration. They are not equilibrated liquid structures and mixed +1/+2 templates are not asserted to be equilibrium solids.

## Authority map

The package now exposes the governing SSDP authority chain explicitly:

| Domain | Canonical owner | Scope |
|---|---|---|
| D1 | `SCIENTIFIC_METHOD.md` | scientific meaning, supported chemistry, normalization, validity, limitations |
| D2 | `NUMERICAL_METHOD.md` | equations, apportionment, finite-cell conversion, numerical/reproducibility semantics |
| D3 | `ARCHITECTURE.md` | ownership, data flow, geometry/serialization boundaries, compatibility architecture |
| D4 | `SPECIFICATION.md` | CLI, defaults, diagnostics, output formats, failure behavior |
| External evidence | `REFERENCES.md` | exact literature/standards locators and source-to-claim mapping |

`LiCl_NaCl_KCl_cell_builder_manual.md` and its PDF are derived technical guides. `tests/test_builder.py` is executable evidence, not scientific authority.

## Scientific contract

For integer cation counts `N_i`, exact chloride count is

```text
N_Cl = sum_i z_i N_i
```

with `z=1` for Group-I cations and `z=2` for Group-II cations. Divalent components are `MCl2`, so their formula mass is `M_M + 2 M_Cl`.

The mixture molar-volume basis is **per mole of cations / salt formula-unit components**, not per total ion count. Pure component volume is

```text
V_i(T) = M_i / rho_i(T)
```

and the ideal mixture baseline is

```text
V_ideal(T,x) = sum_i x_i V_i(T)
```

where `x_i` is the achieved finite-cell cation fraction.

The historical LiCl-NaCl-KCl excess-volume field is deliberately frozen for compatibility with prior calculations. It combines source-backed binary density relations with the package's historical symmetric pair-sum-squared/Kohler-style weighting. This package-defined ternary interpolation is **not** the literal Fellner et al. ternary equation. Fellner data are retained as independent comparison evidence.

Pairs involving newly added Rb, Cs, Be, Mg, Ca, Sr, or Ba have no invented excess-volume term; those interaction channels use ideal-volume interpolation.

## Temperature validity

The runtime warns when an active empirical density relation is evaluated outside its cited source interval. Warnings do not clamp temperature or modify the model.

With `--ideal-only`, the legacy Li/Na/K binary excess terms are inactive, so their binary-range warnings are suppressed while applicable pure-component warnings remain.

The exact source equations and ranges are recorded in `source_density_equations.csv`, `NUMERICAL_METHOD.md`, and the technical manual.

## Quick examples

Legacy LiCl-KCl:

```bash
python build_mixed_alkali_chloride_rocksalt.py \
  --nx 3 --ny 3 --nz 3 \
  --temp-c 800 \
  --x-li 0.585 --x-k 0.415 \
  --format all \
  --output-prefix LiClKCl_eutectic
```

Mixed monovalent/divalent composition:

```bash
python build_mixed_alkali_chloride_rocksalt.py \
  --nx 2 --ny 2 --nz 2 \
  --temp-k 1090 \
  --x-li 0.25 --x-cs 0.25 --x-mg 0.25 --x-ca 0.25 \
  --format all \
  --lammps-atom-style charge \
  --output-prefix LiCsMgCaCl
```

Ideal-only comparison:

```bash
python build_mixed_alkali_chloride_rocksalt.py \
  --nx 2 --ny 2 --nz 2 \
  --temp-k 1000 \
  --x-li 0.50 --x-na 0.25 --x-k 0.25 \
  --ideal-only \
  --format all
```

## Hydrostatic volume strain

`--hydrostatic-strain` (alias `--volume-strain`) is applied after the equilibrium liquid volume:

```text
V = V_eq * (1 + eps_V)
```

A bare value is fractional volume strain (`0.05` = +5% volume); a percent suffix is also accepted (`+5%`). The isotropic length factor is `(1+eps_V)^(1/3)`.

## Outputs

The script writes POSCAR, CIF, and/or LAMMPS data. LAMMPS supports `atomic` or `charge` atom styles; `charge` emits formal `+1`, `+2`, and `-1` charges.

The public serialization contract is defined in `SPECIFICATION.md`; external format references are in `REFERENCES.md`.

## Reproducibility contract

For Li/Na/K compositions with excess terms enabled, the default volume model is frozen to the historical compatibility field. The packaged tests contain binary and ternary numerical anchors and assembled diagnostic oracles. Those anchors protect numerical compatibility; they are not substitutes for empirical scientific validation.

The executable version remains 1.3.4. The SSDP authority/documentation rebuild changes documentation and source comments/provenance routes only; it does not intentionally change executable Python semantics.
