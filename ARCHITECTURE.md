# Architecture Manual - Group-I/Group-II chloride initial melt-cell builder

**Authority domain:** D3 software architecture  
**Authority state:** accepted-current reconstruction for executable version 1.3.4  
**Governing protocol:** SSDP 6.6.0  
**Parents:** `SCIENTIFIC_METHOD.md`, `NUMERICAL_METHOD.md`

## 1. Architectural scope

The package is intentionally a small, single-process command-line program. Its architecture exists to preserve a narrow scientific/numerical pipeline rather than to provide a general simulation framework.

The durable flow is:

```text
CLI inputs
  -> canonical temperature and composition
  -> exact integer cation/Cl counts
  -> thermodynamic molar volume
  -> finite-cell volume and lattice
  -> deterministic randomized site population
  -> POSCAR / CIF / LAMMPS serialization
  -> observational diagnostics
```

No database, service, cache, checkpoint, network dependency, or concurrent worker is part of the runtime architecture.

## 2. Ownership boundaries

| Architectural owner | Responsibility | Must not own |
|---|---|---|
| composition parsing/apportionment | produce one achieved integer cation state | scientific redefinition of component basis |
| thermodynamic functions | concretize D1/D2 pure, ideal, and legacy excess equations | geometry or file-format policy |
| geometry builder | populate the fixed base and auxiliary fractional-site sets | thermodynamic volume fitting |
| lattice builder | map target volume to orthorhombic box lengths | change modeled volume |
| writers | serialize the already constructed state | recompute scientific quantities |
| warning path | report empirical-range status | clamp/refit/change state |
| packaged tests | verify behavior and oracles | act as scientific authority |

The single Python module is a deployment choice, not an assertion that these conceptual ownership boundaries may be collapsed semantically.

## 3. Canonical state and data flow

The achieved integer cation counts are the authoritative runtime composition. Cation fractions are derived once from those counts and reused for thermodynamics and reporting. Chloride count is derived from formal charge, and the generated atom list is checked against that required count.

The thermodynamic calculation precedes and determines the lattice scale. Geometry generation does not feed back into the volume model. File writers consume atoms and lattice vectors but do not own or modify them.

Warnings execute after the scientific state is already determined and therefore remain observational.

## 4. Geometry architecture

The base coordinate template is the conventional rock-salt/halite B1 construction [R6]. Monovalent-only systems use four cation and four chloride sites per conventional cell.

For divalent compositions, the architecture adds one auxiliary chloride candidate per cation site. The number actually occupied is exactly the additional chloride count required by charge neutrality. These auxiliary sites are a package construction mechanism for an initial configuration; they are not a separate claimed equilibrium crystal prototype.

All coordinates are stored as fractional coordinates in the periodic cell until serialization.

## 5. Reproducibility and randomness

A single user-visible integer seed controls both cation-site shuffling and auxiliary-chloride site sampling. Given the same Python/runtime semantics, inputs, and seed, the package is intended to reproduce the same coordinate assignment.

Scientific quantities that precede random assignment - counts, charge, molar volume, density, and lattice dimensions - are seed independent.

## 6. Serialization boundaries

The package supports three file families:

- POSCAR, following the VASP structural input format [R7];
- CIF, using a P1 representation compatible with the CIF standard [R8];
- LAMMPS data, with `atomic` or `charge` atom styles as documented by LAMMPS [R9].

Serialization is downstream of state construction. Writers may control ordering and text formatting but must not change counts, coordinates, charges, masses, or lattice semantics.

## 7. Compatibility architecture

The legacy Li/Na/K public behavior is a protected compatibility boundary:

- historical two-of-three composition completion remains available;
- the frozen Li/Na/K volume field remains the default excess-enabled model;
- legacy formula masses remain fixed;
- monovalent-only geometry remains the original rock-salt template;
- output defaults remain stable unless deliberately versioned.

A future direct Fellner or other alternative thermodynamic model must have a distinct model identity and explicit selection path. It must not replace the legacy model behind the existing default interface.

## 8. Failure behavior

Invalid dimensions, composition conflicts, invalid fractions/counts, nonpositive absolute temperature, invalid hydrostatic strain, or internal charge/count inconsistencies fail synchronously before successful completion. The architecture does not silently repair scientifically invalid inputs with hidden fallbacks.

Empirical range excursions are different: they are valid software inputs but scientifically extrapolative states, so they emit warnings rather than hard failure.

## 9. Non-goals

The architecture does not attempt to provide:

- equilibrium liquid structure generation;
- force-field or molecular-dynamics execution;
- automatic empirical model selection/refitting;
- generalized activity-coefficient or equation-of-state databases;
- persistence/restart or distributed execution.

## 10. D3 -> D4 handoff

`SPECIFICATION.md` owns the stable command-line and serialization behavior. The Python file may refactor internally only while preserving the upstream scientific/numerical contracts and these ownership/data-flow boundaries.

## References

Exact locators are in `REFERENCES.md`. Architecture-level external format/prototype dependencies are [R6]-[R9].
