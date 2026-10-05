# Numerical & Algorithmic Method Paper - Group-I/Group-II chloride initial melt-cell builder

**Authority domain:** D2 numerical/algorithmic semantics  
**Authority state:** accepted-current reconstruction for executable version 1.3.4  
**Governing protocol:** SSDP 6.6.0  
**Parent authority:** `SCIENTIFIC_METHOD.md`

## 1. Inputs and canonical state

The numerical state used for thermodynamics is the achieved integer cation composition of a finite cell, not the unrounded requested fractions.

For positive integers `nx`, `ny`, `nz`, the number of cation sites is

\[
N=4n_xn_yn_z.
\]

A run supplies either cation mole fractions or explicit cation counts, never both.

## 2. Fraction-to-count apportionment

When fractions are supplied, the target real counts are

\[
r_i=Nx_i.
\]

The implementation uses largest-remainder apportionment:

1. `b_i = floor(r_i)`;
2. `R = N - sum_i b_i`;
3. add one count to the `R` species with largest fractional remainder `r_i-b_i`;
4. exact ties are resolved by the fixed cation order `Li, Na, K, Rb, Cs, Be, Mg, Ca, Sr, Ba`.

The achieved fractions are then

\[
\hat{x}_i=\frac{N_i}{N}.
\]

All thermodynamic mixture calculations use `hat{x}_i` so composition, mass, charge, cell volume, and serialized species counts refer to the same finite state.

The legacy Li/Na/K parser retains its historical completion rule: if exactly two of the three Li/Na/K fractions (or counts) are supplied, the third is completed to unity (or to `N`). When any extended cation is used, omitted species are zero and supplied fractions must sum to one.

## 3. Formula masses and pure volumes

For monovalent `M+`,

\[
M_{MCl}=M_M+M_{Cl}.
\]

For divalent `M^{2+}`,

\[
M_{MCl_2}=M_M+2M_{Cl}.
\]

The legacy LiCl, NaCl, and KCl formula masses are frozen compatibility constants. Other packaged formula masses are calculated from the listed atomic masses. `REFERENCES.md` [R11] provides formula/molecular-weight provenance, but frozen legacy constants are not re-normalized to newer tables.

Each pure molar volume is

\[
V_i^0(T)=\frac{M_i}{\rho_i(T)}.
\]

## 4. Pure-density equations

### 4.1 Legacy Li/Na/K relations

The runtime relations are

\[
\rho_{LiCl}=1.8721-4.243\times10^{-4}T_K,
\]

\[
\rho_{NaCl}=2.1390-5.444\times10^{-4}T_K,
\]

and pure KCl is the NaCl-KCl fitted surface evaluated at `C=100 mol% KCl`:

\[
\rho=2.14250-5.45436\times10^{-4}T_K
+2.79726\times10^{-8}C^3
-9.89901\times10^{-7}T_KC
+2.78444\times10^{-10}CT_K^2.
\]

The KCl-LiCl fitted surface is

\[
\rho=1.98715-5.89069\times10^{-4}T_K
+2.02443\times10^{-3}C
+5.31098\times10^{-11}T_K^3
+1.21294\times10^{-7}C^3
-2.16826\times10^{-8}T_KC^2,
\]

with `C=100 y_LiCl`. These relations and their source interpretation are supported by [R1] and [R2].

### 4.2 Additional pure salts

For RbCl, CsCl, BeCl2, MgCl2, CaCl2, SrCl2, and BaCl2, the Janz 1988 representation [R12] is

\[
\rho(T_C)=\rho_m-k(T_C-T_m).
\]

The exact parameter table used by the runtime is in `source_density_equations.csv` and is reproduced in the technical manual.

## 5. Ideal mixture volume

On the D1 component basis,

\[
V_m^{ideal}(T,\hat{\mathbf{x}})=\sum_i \hat{x}_iV_i^0(T).
\]

This expression is the complete model when `--ideal-only` is selected.

## 6. Historical Li/Na/K binary excess construction

For an available binary pair `i-j`, define

\[
s_{ij}=\hat{x}_i+\hat{x}_j,
\qquad
y_i=\frac{\hat{x}_i}{s_{ij}},
\qquad
y_j=\frac{\hat{x}_j}{s_{ij}}.
\]

The source binary density gives

\[
V_{ij}^{bin}=\frac{y_iM_i+y_jM_j}{\rho_{ij}(T,y)},
\]

and the binary excess volume is

\[
V_{ij}^{E,bin}=V_{ij}^{bin}-y_iV_i^0-y_jV_j^0.
\]

For LiCl-NaCl, the historical implementation reconstructs the arbitrary-composition binary from the equimolar source relation. With

\[
L_0=4\left(V_{50:50}-V_{50:50}^{ideal}\right),
\]

it uses

\[
V_{LiNa}^{bin}=V_{LiNa}^{ideal}+y_{Li}y_{Na}L_0.
\]

This is a package-specific single-anchor binary construction; [R5] is general excess-property provenance, not a claim of a fitted higher-order Redlich-Kister series.

## 7. Historical ternary/general-mixture assembly

Let

\[
\mathcal{K}=\{LiCl-NaCl,\ NaCl-KCl,\ KCl-LiCl\}.
\]

With excess terms enabled, the accepted numerical model is

\[
V_m(T,\hat{\mathbf{x}})=
\sum_i\hat{x}_iV_i^0(T)
+\sum_{(i,j)\in\mathcal{K}}
(\hat{x}_i+\hat{x}_j)^2
V_{ij}^{E,bin}(T,y_i,y_j).
\]

The pair-sum-squared weighting is the historical symmetric Kohler-style construction [R4]. This exact equation, not the name "Kohler" by itself, is the normative compatibility target.

No excess term is added for an unsupported pair.

## 8. Mixture mass, density, and finite-cell conversion

The mixture formula-component molar mass is

\[
M_{mix}=\sum_i\hat{x}_iM_i,
\]

and the reported equilibrium density is

\[
\rho_{mix}=\frac{M_{mix}}{V_m}.
\]

With exact SI Avogadro constant `N_A = 6.02214076e23 mol^-1` [R10],

\[
V_{cell}^{eq}=N\frac{V_m}{N_A}10^{24}\ \mathrm{\mathring{A}}^3.
\]

For `n_conv = nx ny nz`, the equivalent conventional-cell edge used for the orthorhombic supercell is

\[
a_{conv}=\left(\frac{V_{cell}}{n_{conv}}\right)^{1/3}.
\]

The lattice lengths are `nx*a_conv`, `ny*a_conv`, and `nz*a_conv`, whose determinant equals the target total volume up to ordinary binary64 rounding and decimal serialization.

## 9. Hydrostatic volume strain

After the equilibrium thermodynamic calculation,

\[
V_m^{strained}=V_m^{eq}(1+\epsilon_V),
\]

with required `epsilon_V > -1`. The isotropic length factor is

\[
s=(1+\epsilon_V)^{1/3}.
\]

This is a volume-strain convention. The strain changes the generated box but does not refit or redefine the equilibrium mixture model.

## 10. Temperature-range diagnostics

Source validity intervals are inclusive. The comparison predicate uses an absolute tolerance of `1e-9` degree solely to prevent false endpoint warnings caused by binary floating-point conversion between Kelvin and Celsius.

Legacy pure correlations are diagnosed whenever their species is active. Legacy binary diagnostics are evaluated only when excess terms are enabled and both pair members are materially active. `--ideal-only` therefore suppresses binary-model validity warnings while preserving active pure-component warnings.

The warning system is observational: it has no data-flow path that changes the temperature, composition, molar volume, cell volume, or coordinates.

## 11. Randomized site assignment

After counts and volume are fixed, species are assigned to the fixed coordinate template with Python's deterministic pseudo-random generator initialized from the public integer seed. The default seed is `12345`.

Randomization affects which species occupies which template site and which auxiliary chloride sites are selected; it does **not** affect species counts, charge, thermodynamic volume, or lattice dimensions.

The exact pseudo-random implementation is a D4 reproducibility mechanism, not a D1 scientific model.

## 12. Numerical precision and equivalence

The implementation uses Python binary64 floating-point arithmetic. Discrete invariants - species counts, chloride count, and formal charge construction - are exact integer identities. Scientific floating observables are expected to agree with independent formulas to binary64/printed-roundoff precision.

The legacy Li/Na/K field is a compatibility invariant: accidental numerical changes to its equations, coefficients, normalization, or pair weighting are not permitted without an explicitly different model identity.

## 13. D2 -> D3 handoff

D3/D4 must preserve:

- a single canonical achieved integer composition for thermodynamics and serialization;
- exact integer charge construction;
- the D1/D2 component normalization;
- the stated pure/binary correlations and frozen legacy assembly equation;
- deterministic finite-cell volume conversion;
- hydrostatic strain as a post-equilibrium scale operation;
- diagnostics that remain observational;
- explicit separation of compatibility oracles from empirical-validity evidence.

## References

Exact source identities and locators are in `REFERENCES.md`, especially [R1], [R2], [R4], [R5], [R10], [R11], and [R12].
