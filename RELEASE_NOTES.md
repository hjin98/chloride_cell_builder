# SSDP authority/documentation rebuild - executable version 1.3.4

This distribution rebuilds the human-facing documentation around an explicit SSDP 6.6.0 authority chain while retaining the reviewed version-1.3.4 executable semantics.

- adds `SCIENTIFIC_METHOD.md` (D1), `NUMERICAL_METHOD.md` (D2), `ARCHITECTURE.md` (D3), `SPECIFICATION.md` (D4), and `REFERENCES.md`;
- rewrites the README and technical manual as current-system documents that route to those canonical owners instead of using release-history prose as authority;
- expands exact source locators for scientific data, SI constants, crystallographic provenance, and serialization formats;
- updates only comments/provenance routing in the Python source; the parsed Python AST is required to remain identical to the reviewed v1.3.4 implementation;
- regenerates the PDF manual, manifest, checksums, and distributable archive.

No thermodynamic equation, coefficient, normalization rule, charge rule, integer apportionment rule, geometry rule, CLI default, warning predicate, or test oracle is intentionally changed by this documentation rebuild.

---

# Referenced release notes - version 1.3.4

This patch release strengthens the packaged regression oracle for the version-1.3.3 validity-diagnostic repair without changing thermodynamic, stoichiometric, volume, geometry, or runtime diagnostic behavior.

Key changes:

- **no thermodynamic volume equation, coefficient, composition rule, charge rule, warning predicate, box-scaling rule, or coordinate construction was changed** relative to version 1.3.3;
- adds an assembled-CLI positive oracle for all three legacy LiCl-NaCl, NaCl-KCl, and KCl-LiCl binary validity warnings with the excess model enabled;
- pairs that positive oracle with the existing assembled `--ideal-only` negative oracle, so the suite now discriminates both under-warning and over-suppression failures in `main()` state propagation;
- mutation-checking the assembled warning call by forcing `use_excess=False` now causes the new regression test to fail, whereas the version-1.3.3 suite could remain green under that defect;
- preserves the frozen Li/Na/K compatibility field, pure-component diagnostics, endpoint behavior, mixed-valence normalization, finite-cell construction, and generated structures unchanged.

Compatibility policy: this is a test-oracle/release-evidence repair only. The executable scientific and diagnostic contract is unchanged from version 1.3.3 except for the reported package version.

---

# Referenced release notes - version 1.3.3

This patch release repairs the final validity-diagnostic mismatch found in the independent v1.3.2 closure review while preserving all thermodynamic and structural numerical behavior.

Key changes:

- **no thermodynamic volume equation, coefficient, composition rule, charge rule, box-scaling rule, or coordinate construction was changed** relative to version 1.3.2;
- makes legacy LiCl-NaCl, NaCl-KCl, and KCl-LiCl binary-range warnings conditional on those binary excess terms actually being enabled;
- under `--ideal-only`, retains warnings for active pure-component density correlations while suppressing validity warnings for inactive legacy binary excess models;
- adds a direct diagnostic regression test and an assembled CLI negative oracle for the `--ideal-only` path;
- preserves the frozen Li/Na/K compatibility field and all version-1.3.2 empirical endpoint behavior unchanged.

Compatibility policy: this is an observational diagnostic repair only. It does not alter any requested thermodynamic state, calculated molar volume, finite-cell volume, density, box dimensions, species counts, or generated coordinates.

---

# Referenced release notes - version 1.3.2

This patch release repairs the remaining empirical-validity and regression-oracle issues identified by the independent v1.3.0/v1.3.1 numerical review while preserving the frozen LiCl-NaCl-KCl volume field exactly.

Key changes:

- **no thermodynamic volume equation or coefficient was changed** relative to version 1.3.1;
- adds source-range warnings for every active LiCl, NaCl, and KCl pure-component baseline correlation, including when it is used inside a mixture;
- uses the reported 893-1123 K range of the 50:50 LiCl-NaCl source correlation for the legacy LiCl-NaCl excess-term diagnostic;
- records broad source-study envelopes for the NaCl-KCl and KCl-LiCl fitted surfaces and labels their exact liquid-data interval as composition dependent;
- fixes false out-of-range warnings at exact Celsius endpoints caused only by floating-point K-to-C conversion roundoff;
- strengthens tests with independent hard-coded source anchors rather than deriving expected densities from production coefficient tables;
- independently reconstructs an ideal mixed-valence molar-volume case instead of using production pure-volume helpers as its oracle;
- decouples the LAMMPS charge-serialization test from the production thermodynamic volume path;
- adds explicit endpoint/just-outside checks for all ten pure correlations and all three legacy binary diagnostics;
- adds an end-to-end CLI finite-cell-volume anchor for the historical 3x3x3 LiCl-KCl case.

Compatibility policy: temperature diagnostics are observational only. They do not clamp, reject, refit, or otherwise change the requested temperature, composition, volume, box dimensions, or generated coordinates. The version-1.3.1 Li/Na/K reproducibility contract remains in force.

---

# Referenced release notes - version 1.3.1

This patch release reconciles the provenance and reproducibility contract of the inherited LiCl-NaCl-KCl volume model without changing its numerical behavior.

Key changes:

- **no thermodynamic volume equation or coefficient was changed** relative to version 1.3.0;
- designates the existing LiCl-NaCl-KCl pairwise excess-volume plus symmetric Kohler-style construction as the **legacy Li/Na/K compatibility model**;
- clarifies that the pure and binary-edge density relations are source-backed, but the assembled ternary-interior formula is the package's historical interpolation and is **not** the direct Fellner et al. (1980) ternary excess-volume equation;
- retains Fellner et al. (1980) as external validation/comparison evidence rather than presenting it as the equation implemented by the default model;
- documents the existing Fellner-point discrepancies transparently instead of refitting the model, because prior scientific work may depend on the historical volume field;
- adds multiple frozen ternary-interior regression anchors so future maintenance cannot silently alter the legacy Li/Na/K volume behavior;
- clarifies that those frozen anchors are compatibility oracles, not empirical-validity oracles.

Compatibility policy: any future implementation of a direct Fellner or other alternative ternary model should be separately named and explicitly selected. It must not silently replace the legacy default because doing so would break numerical reproducibility of prior Li/Na/K calculations.

---

# Referenced release notes - version 1.3.0

This release generalizes the LiCl-NaCl-KCl builder to stable Group-I and Group-II metal chlorides with source-backed pure-melt density models.

Key changes:

- added `Rb+`, `Cs+`, `Be2+`, `Mg2+`, `Ca2+`, `Sr2+`, and `Ba2+` while retaining `Li+`, `Na+`, and `K+`;
- generalized composition parsing, largest-remainder finite-cell rounding, species ordering, masses, formal charges, and output writers;
- enforces exact charge balance with `N_Cl = sum_i z_i N_i`, so mixed +1/+2 systems contain the correct number of chloride ions;
- preserves the original rocksalt site topology for monovalent-only cells and adds charge-required chloride ions on an auxiliary tetrahedral-like site set when divalent cations are present;
- generalized the volume model to interpolate pure molten `ACl`/`BCl2` molar volumes at the achieved cation composition and requested temperature;
- retained the original source-backed LiCl-NaCl-KCl binary excess-volume/Kohler model exactly on its original composition subspace;
- leaves newly introduced binary interactions ideal rather than inventing unsupported excess-volume parameters;
- added Janz 1988 pure-melt density correlations and validity-range warnings for RbCl, CsCl, BeCl2, MgCl2, CaCl2, SrCl2, and BaCl2;
- preserved the legacy Li/Na/K density equations and legacy Li/Na/K CLI completion behavior for backward compatibility;
- added focused regression tests for legacy Li/Na/K values, pure-density anchors, composition interpolation, +1/+2 charge balance, and generated LAMMPS charges.

Scope note: HCl is not treated as the same simple metal-chloride ionic melt, and this package does not fabricate density models for the radioactive FrCl or RaCl2 endpoints. Therefore the supported set is the stable alkali/alkaline-earth metal cations with compatible source data.

---

# Referenced release notes - version 1.2.0

This release adds optional hydrostatic strain on top of the temperature- and composition-dependent equilibrium liquid volume.

Key changes:

- added `--hydrostatic-strain` with alias `--volume-strain`;
- defined the option using a **volume-strain convention**, `V=V_eq(1+eps_V)`;
- accepts either fractional or percent-suffixed input;
- supports positive expansion and negative compression while rejecting strains at or below `-100%`;
- applies the exact isotropic lattice factor `(1+eps_V)^(1/3)`;
- preserves legacy geometry when the option is omitted.

---

# Referenced release notes - version 1.1.0

This release consolidated the thermodynamic theory, implementation details, source equations, validation examples, and structural-file examples into one package and documented the direct binary-density-derived excess-volume construction and Kohler ternary extrapolation.
