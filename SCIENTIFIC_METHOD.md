# Scientific Method Paper - Group-I/Group-II chloride initial melt-cell builder

**Authority domain:** D1 scientific/mathematical meaning  
**Authority state:** accepted-current reconstruction for executable version 1.3.4  
**Governing protocol:** SSDP 6.6.0  
**Change status:** no new scientific model is introduced here. This document reconstructs the already accepted contract supplied for the v1.3.4 review and separates it from numerical and software detail.

## 1. Scientific purpose

The package constructs a periodic **initial configuration** for subsequent simulation of molten chloride mixtures. Its scientific role is to provide:

1. an exactly charge-neutral ionic composition;
2. a finite-cell volume consistent with the package's accepted temperature/composition molar-volume model; and
3. a compact periodic coordinate template suitable for later melting/equilibration.

The generated coordinates are **not** an equilibrated liquid structure and, for mixed monovalent/divalent systems, are **not** asserted to be an equilibrium solid phase.

## 2. Supported chemical domain

The supported practical cations are

- monovalent: Li+, Na+, K+, Rb+, Cs+;
- divalent: Be2+, Mg2+, Ca2+, Sr2+, Ba2+.

Each monovalent component is interpreted as `MCl`; each divalent component is interpreted as `MCl2`. H+, Fr+, and Ra2+ are outside the supported volume-model domain: HCl is not treated as the same simple metal-chloride ionic melt, and the packaged reference set does not provide compatible critically evaluated pure-melt density correlations for FrCl or RaCl2.

## 3. Composition and charge invariant

For integer cation counts `N_i` and formal cation valence `z_i in {1,2}`, the chloride count is

\[
N_{Cl}=\sum_i z_i N_i.
\]

The formal ionic charge therefore satisfies

\[
\frac{Q}{e}=\sum_i z_iN_i-N_{Cl}=0
\]

exactly. Charge neutrality is a scientific invariant, not a post-hoc numerical tolerance.

## 4. Thermodynamic observable and normalization

The primary modeled observable is the equilibrium mixture molar volume

\[
V_m(T,\mathbf{x})
\]

in `cm^3/mol-cation`, where `x_i` is the achieved **cation mole fraction**. Equivalently, because every `MCl` or `MCl2` formula unit contains one cation, this is a per-mole-of-salt-formula-components basis.

The basis is **not** per mole of all ions. This distinction is essential for mixtures containing `MCl2`.

For each pure salt component `i`, the pure molten molar volume is

\[
V_i^0(T)=\frac{M_i}{\rho_i^0(T)},
\]

where `M_i` is the formula mass of `MCl` or `MCl2` and `rho_i^0` is the applicable pure-melt density correlation.

## 5. Scientific mixture model

### 5.1 Ideal baseline

All supported mixtures have the ideal-volume baseline

\[
V_m^{ideal}(T,\mathbf{x})=\sum_i x_iV_i^0(T).
\]

For interactions involving newly added Rb, Cs, Be, Mg, Ca, Sr, or Ba components, no unsourced excess-volume term is invented. Such interaction channels are therefore treated by ideal-volume interpolation in the current model.

This is an **approximation/closure assumption**, not a scientific assertion that the real excess volume is zero.

### 5.2 Frozen LiCl-NaCl-KCl compatibility model

For the historical LiCl-NaCl-KCl subspace, the package deliberately preserves the established pairwise excess-volume field used by earlier versions and earlier scientific calculations. Source-backed binary density information from [R1] and [R2] is converted to binary excess molar volume and assembled using the historical symmetric Kohler-style weighting described normatively in `NUMERICAL_METHOD.md`.

The resulting ternary-interior field is a **package-defined compatibility interpolation**. It is not the literal ternary equation reported by Fellner et al. [R3].

Preserving the historical Li/Na/K numerical field is an accepted compatibility requirement. It must not be silently replaced merely because another direct ternary model may be scientifically preferable for a different purpose.

### 5.3 `--ideal-only`

`--ideal-only` removes the historical Li/Na/K excess contributions and evaluates only the pure-component ideal-volume baseline. It does not change composition, charge balance, pure-component density evaluation, or the requested temperature.

## 6. Empirical sources and validity

The legacy Li/Na/K density relations and source envelopes are supported by Janz et al. [R1] and Van Artsdalen-Yaffe [R2]. The seven additional pure chloride correlations are taken from Janz's critically evaluated 1988 compilation [R12]. Exact equations and ranges are recorded in `source_density_equations.csv` and `NUMERICAL_METHOD.md`.

Evaluation outside a cited source interval is permitted for compatibility and continuity, but is scientifically interpreted as **extrapolation**. Runtime warnings are observational diagnostics: they do not clamp, reject, or refit the requested thermodynamic state.

For Li/Na/K ternary mixtures, the Fellner data [R3] provide an external comparison set. They are retained specifically to expose the empirical discrepancy of the frozen compatibility field rather than to redefine that field.

## 7. Finite-cell interpretation

For `N_cation` cations, the equilibrium periodic cell volume is derived from the mixture molar volume as

\[
V_{cell}^{eq}=N_{cation}\frac{V_m}{N_A}10^{24}\;\mathrm{\mathring{A}}^3,
\]

with the exact SI Avogadro constant [R10].

An optional hydrostatic volume strain may scale this equilibrium volume after the thermodynamic prediction. The strain does not alter the underlying equilibrium volume model.

## 8. Scientific limitations and uncertainty

The principal model limitations are:

- the frozen Li/Na/K ternary interior is a compatibility interpolation, not the direct Fellner ternary equation;
- no comprehensive excess-volume database is supplied for all newly added binary or higher-order interactions;
- source correlations have finite empirical temperature ranges, so out-of-range use is extrapolative;
- the generated geometry is an initial construction template, not an equilibrium structural model.

Accordingly, the package is suitable for creating volume-scaled starting cells. It is not presented as a universal high-accuracy multicomponent chloride equation of state.

## 9. External adequacy and falsification

Scientific adequacy is checked by separating three questions:

1. **Source fidelity:** the packaged pure/binary correlations and ranges must agree with [R1], [R2], and [R12].
2. **Compatibility fidelity:** the Li/Na/K numerical field must remain identical to its accepted historical implementation unless a separately named model is introduced.
3. **Empirical comparison:** independent ternary data such as [R3] may quantify discrepancy but do not silently redefine the compatibility default.

A future claim of improved physical accuracy requires a separately identified scientific model, its own validity basis, and explicit selection; it must not masquerade as maintenance of the existing compatibility model.

## 10. D1 -> D2 handoff

D2 must preserve the following invariants:

- component basis is per mole of cations / salt formula units;
- `MCl2` formula masses contain two chlorides;
- exact formal charge neutrality is `N_Cl = sum_i z_i N_i`;
- pure volumes are `M_i/rho_i(T)` with the cited correlation identities;
- unsupported new-species interaction channels remain ideal in the current model;
- the legacy Li/Na/K field remains numerically frozen;
- source-range diagnostics remain observational and do not alter the requested state;
- finite-cell volume follows the stated Avogadro conversion.

## References

Reference labels are defined with exact locators in `REFERENCES.md`. The principal D1 sources are [R1], [R2], [R3], [R10], [R11], and [R12].
