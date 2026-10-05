---
title: "Group-I / Group-II Chloride Initial Melt-Cell Builder"
subtitle: "Technical Manual - executable version 1.3.4"
author: "SSDP authority/documentation reconstruction"
date: "2026-10-04"
geometry: margin=0.8in
fontsize: 10pt
header-includes:
  - \usepackage{booktabs}
  - \usepackage{longtable}
  - \usepackage{microtype}
---

# 1. Purpose, status, and authority

`build_mixed_alkali_chloride_rocksalt.py` constructs periodic starting cells for molten-chloride simulations. The output is designed for subsequent melting/equilibration; it is not an equilibrated liquid structure and the mixed +1/+2 coordinate template is not an asserted equilibrium solid.

This manual is a **derived technical guide**. The canonical SSDP authority set is:

- D1 scientific meaning: `SCIENTIFIC_METHOD.md`;
- D2 numerical method: `NUMERICAL_METHOD.md`;
- D3 architecture: `ARCHITECTURE.md`;
- D4 public contract: `SPECIFICATION.md`;
- exact external source locators: `REFERENCES.md`.

The executable remains version 1.3.4. This documentation rebuild does not introduce a new thermodynamic, stoichiometric, geometry, or compatibility model.

# 2. Supported chemical domain

Supported cations are:

- Group I: Li+, Na+, K+, Rb+, Cs+;
- Group II: Be2+, Mg2+, Ca2+, Sr2+, Ba2+.

Monovalent components are `MCl`; divalent components are `MCl2`. H+, Fr+, and Ra2+ are outside the packaged volume-model domain for the reasons stated in `SCIENTIFIC_METHOD.md`.

# 3. Finite composition and exact charge balance

For `nx`, `ny`, and `nz` conventional template cells, the number of cation sites is

$$
N=4n_xn_yn_z.
$$

The CLI accepts either cation mole fractions or explicit cation counts. Fraction inputs are mapped to integer counts by largest-remainder apportionment; all thermodynamic calculations then use the achieved finite-cell fractions

$$
\hat{x}_i=N_i/N.
$$

For formal cation valence $z_i\in\{1,2\}$, chloride count is

$$
N_{Cl}=\sum_i z_iN_i,
$$

so

$$
Q/e=\sum_i z_iN_i-N_{Cl}=0
$$

exactly.

The historical Li/Na/K parser retains its two-of-three completion behavior. When any extended cation is used, omitted species are zero and the supplied fractions must sum to one (or supplied counts must sum to $N$).

# 4. Formula masses and normalization

For a monovalent component,

$$
M_{MCl}=M_M+M_{Cl},
$$

and for a divalent component,

$$
M_{MCl_2}=M_M+2M_{Cl}.
$$

The thermodynamic basis is **per mole of cations / per mole of salt formula-unit components**. It is not per mole of all ions.

Legacy LiCl, NaCl, and KCl numerical masses are deliberately frozen for compatibility with earlier calculations. Other formula masses are constructed from the packaged atomic-mass constants. NIST Chemistry WebBook SRD 69 is retained as a formula/molecular-weight reference [11], but the frozen legacy numbers are compatibility parameters rather than a claim to be the newest recommended atomic-weight table.

# 5. Temperature convention

The CLI accepts either Kelvin or degrees Celsius:

$$
T_K=T_C+273.15.
$$

Absolute temperature must be positive. Empirical source-range checks are inclusive. A numerical tolerance of $10^{-9}$ degree is applied only to the range predicate to avoid false warnings at exact K/C endpoints; it does not change model temperature or widen the physical source range.

# 6. Pure-component molten volumes

For every supported salt component $i$,

$$
V_i^0(T)=\frac{M_i}{\rho_i(T)}.
$$

## 6.1 Legacy LiCl, NaCl, and KCl relations

From the legacy Janz/NIST compilation [1], the runtime uses

$$
\rho_{LiCl}=1.8721-4.243\times10^{-4}T_K,
$$

$$
\rho_{NaCl}=2.1390-5.444\times10^{-4}T_K.
$$

The NaCl-KCl fitted surface is

$$
\begin{aligned}
\rho={}&2.14250-5.45436\times10^{-4}T_K
+2.79726\times10^{-8}C^3\\
&-9.89901\times10^{-7}T_KC
+2.78444\times10^{-10}CT_K^2,
\end{aligned}
$$

with $C=100y_{KCl}$. Pure KCl is evaluated at $C=100$.

The KCl-LiCl fitted surface [1,2] is

$$
\begin{aligned}
\rho={}&1.98715-5.89069\times10^{-4}T_K
+2.02443\times10^{-3}C\\
&+5.31098\times10^{-11}T_K^3
+1.21294\times10^{-7}C^3
-2.16826\times10^{-8}T_KC^2,
\end{aligned}
$$

with $C=100y_{LiCl}$.

The LiCl-NaCl equimolar source relation is

$$
\rho_{50:50}=2.0121-4.858\times10^{-4}T_K.
$$

## 6.2 Added pure chloride correlations

For RbCl, CsCl, BeCl2, MgCl2, CaCl2, SrCl2, and BaCl2, the critically evaluated Janz 1988 form [12] is

$$
\rho(T_C)=\rho_m-k(T_C-T_m).
$$

| Salt | $T_m$ (deg C) | $\rho_m$ (g cm$^{-3}$) | $k$ (g cm$^{-3}$ deg C$^{-1}$) | Upper limit (deg C) |
|---|---:|---:|---:|---:|
| RbCl | 715 | 2.248 | 0.000883 | 923 |
| CsCl | 645 | 2.790 | 0.001065 | 906 |
| BeCl2 | 415 | 1.540 | 0.001100 | 473 |
| MgCl2 | 714 | 1.680 | 0.000271 | 826 |
| CaCl2 | 775 | 2.085 | 0.000422 | 950 |
| SrCl2 | 874 | 2.727 | 0.000578 | 1037 |
| BaCl2 | 962 | 3.174 | 0.000681 | 1081 |

The machine-readable source table is `source_density_equations.csv`.

# 7. Ideal mixture baseline

Using achieved cation fractions,

$$
V_m^{ideal}(T,\hat{\mathbf{x}})=\sum_i\hat{x}_iV_i^0(T).
$$

This is the complete thermodynamic volume model under `--ideal-only`.

For newly introduced interaction channels involving Rb, Cs, Be, Mg, Ca, Sr, or Ba, no unsourced excess parameter is invented. Ideal interpolation is therefore an explicit approximation, not a claim that the real excess volume is zero.

# 8. Historical Li/Na/K excess-volume compatibility model

For a legacy binary pair $i-j$, define

$$
s_{ij}=\hat{x}_i+\hat{x}_j,
\qquad
y_i=\frac{\hat{x}_i}{s_{ij}},
\qquad
y_j=\frac{\hat{x}_j}{s_{ij}}.
$$

The source binary density gives

$$
V_{ij}^{bin}=\frac{y_iM_i+y_jM_j}{\rho_{ij}(T,y)},
$$

and

$$
V_{ij}^{E,bin}=V_{ij}^{bin}-y_iV_i^0-y_jV_j^0.
$$

For LiCl-NaCl, the historical implementation reconstructs arbitrary composition from the equimolar source relation with

$$
L_0=4\left(V_{50:50}-V_{50:50}^{ideal}\right)
$$

and

$$
V_{LiNa}^{bin}=V_{LiNa}^{ideal}+y_{Li}y_{Na}L_0.
$$

Reference [5] provides general excess-property representation provenance; the exact equation above is the package contract.

The global mixture model with excess enabled is

$$
V_m(T,\hat{\mathbf{x}})=
\sum_i\hat{x}_iV_i^0(T)
+\sum_{(i,j)\in\mathcal{K}}
(\hat{x}_i+\hat{x}_j)^2V_{ij}^{E,bin}(T,y_i,y_j),
$$

where

$$
\mathcal{K}=\{LiCl-NaCl,\ NaCl-KCl,\ KCl-LiCl\}.
$$

The pair-sum-squared construction is the historical symmetric Kohler-style extrapolation [4]. The exact package equation above, rather than the generic method name, defines compatibility.

## 8.1 Fellner comparison and model identity

Fellner, Votava, and Chrenkova-Paucirova [3] measured the LiCl-NaCl-KCl ternary system and reported calculated/experimental comparisons. Their data remain external validation evidence.

The package's ternary field is **not** the direct Fellner ternary equation. It is the historical pairwise/Kohler compatibility interpolation. This distinction is intentional because earlier scientific calculations may depend on its exact numerical values.

At 800 deg C, representative comparison points are:

| $x_{LiCl}$ | $x_{NaCl}$ | $x_{KCl}$ | Legacy model $\rho$ | Fellner experimental $\rho$ | Difference |
|---:|---:|---:|---:|---:|---:|
| 1/6 | 1/6 | 2/3 | 1.469016 | 1.498000 | -1.935% |
| 1/6 | 2/3 | 1/6 | 1.525644 | 1.519000 | +0.437% |
| 1/3 | 1/3 | 1/3 | 1.498322 | 1.492000 | +0.424% |
| 2/3 | 1/6 | 1/6 | 1.497401 | 1.455000 | +2.914% |

These discrepancies are retained as model limitation evidence, not used to refit the legacy default.

# 9. Validity diagnostics

The explicitly tracked legacy ranges are:

| Runtime correlation | Diagnostic interval | Meaning |
|---|---:|---|
| pure LiCl | 880-1070 K | Janz/NIST pure-LiCl tabulation |
| pure NaCl | 1090-1220 K | Janz/NIST pure-NaCl tabulation |
| pure KCl via NaCl-KCl surface | 1060-1200 K | tabulated KCl endpoint interval |
| LiCl-NaCl 50:50 anchor | 893-1123 K | source range used to construct legacy excess term |
| NaCl-KCl fitted surface | 943-1300 K | broad source-study envelope; exact liquid-data interval composition dependent |
| KCl-LiCl fitted surface | 668-1212 K | broad source-study envelope; exact liquid-data interval composition dependent |

Every active pure correlation is checked because it contributes to the ideal baseline. Legacy binary warnings are checked only when the excess model is enabled and both species of that pair are active.

Thus `--ideal-only` suppresses legacy binary-model warnings but not applicable pure-component warnings. A warning indicates extrapolation; it does not alter the requested state or the computed volume.

# 10. Finite-cell volume and hydrostatic strain

The mixture formula-component molar mass is

$$
M_{mix}=\sum_i\hat{x}_iM_i,
$$

and the equilibrium density is

$$
\rho_{mix}=\frac{M_{mix}}{V_m}.
$$

For $N$ cations,

$$
V_{cell}^{eq}=N\frac{V_m}{N_A}10^{24}\ \mathrm{\mathring{A}}^3,
$$

where the SI defines $N_A=6.02214076\times10^{23}\ \mathrm{mol^{-1}}$ exactly [10].

Optional hydrostatic **volume strain** is applied afterward:

$$
V=V_{eq}(1+\epsilon_V),
$$

with isotropic length factor

$$
s=(1+\epsilon_V)^{1/3}.
$$

`epsilon_V` must exceed `-1`.

# 11. Coordinate construction

The base template is the conventional rock-salt/halite B1 structure [6], with four cation and four chloride positions per conventional cell. Monovalent-only compositions retain this topology.

Each divalent cation requires one additional chloride for formal charge balance. The package therefore defines four auxiliary chloride candidate positions per conventional cell at quarter-coordinate tetrahedral-like sites. A seeded random sample of those sites is occupied until the exact required chloride count is reached.

These auxiliary positions are a construction device for an initial cell, not an asserted equilibrium structure.

# 12. Randomization and reproducibility

The public `--seed` defaults to `12345` and controls:

1. shuffling cation species over cation template sites; and
2. selection of auxiliary chloride sites when divalent cations are present.

Species counts, formal charge, thermodynamic volume, and lattice dimensions are fixed before randomization and are therefore seed independent.

# 13. Output formats

The package writes:

- POSCAR using the VASP structural input convention [7];
- CIF using a P1 crystallographic information file representation [8];
- LAMMPS data using `atomic` or `charge` atom styles [9].

With `--lammps-atom-style charge`, per-atom formal charges are written explicitly. Users must select the matching LAMMPS `atom_style` before `read_data` [9].

# 14. Command-line examples

Legacy LiCl-KCl:

```bash
python build_mixed_alkali_chloride_rocksalt.py \
  --nx 3 --ny 3 --nz 3 \
  --temp-c 800 \
  --x-li 0.585 --x-k 0.415 \
  --format all \
  --output-prefix LiClKCl_eutectic
```

Mixed +1/+2:

```bash
python build_mixed_alkali_chloride_rocksalt.py \
  --nx 2 --ny 2 --nz 2 \
  --temp-k 1090 \
  --x-li 0.25 --x-cs 0.25 --x-mg 0.25 --x-ca 0.25 \
  --format all \
  --lammps-atom-style charge \
  --output-prefix LiCsMgCaCl
```

Ideal-only A/B comparison:

```bash
python build_mixed_alkali_chloride_rocksalt.py \
  --nx 2 --ny 2 --nz 2 \
  --temp-k 1000 \
  --x-li 0.50 --x-na 0.25 --x-k 0.25 \
  --ideal-only \
  --format all
```

# 15. Validation and evidence interpretation

The packaged tests verify discrete invariants, density anchors, mixed-valence normalization, finite-cell conversion, output serialization, source-range boundaries, and compatibility fixtures.

The v1.3.4 suite also contains complementary assembled CLI warning oracles:

- active extrapolated legacy binary model -> binary warning must appear;
- `--ideal-only` -> inactive legacy binary warning must not appear while pure warnings remain.

Tests are evidence of the implementation contract. They do not replace the literature sources or the D1/D2 authority.

The dominant remaining scientific limitation is the lack of a general excess-volume model for every newly supported multicomponent interaction. The package should therefore be interpreted as a physically referenced starting-volume builder, not as a universal high-accuracy equation of state.

# References

\small

[1] G. J. Janz, R. P. T. Tomkins, C. B. Allen, J. R. Downey, Jr., G. L. Gardner, U. Krebs, and S. K. Singer, "Molten Salts: Volume 4, Part 2, Chlorides and Mixtures - Electrical Conductance, Density, Viscosity, and Surface Tension Data," *J. Phys. Chem. Ref. Data* **4**(4), 871-1178 (1975). [DOI](https://doi.org/10.1063/1.555527)

[2] E. R. Van Artsdalen and I. S. Yaffe, "Electrical Conductance and Density of Molten Salt Systems: KCl-LiCl, KCl-NaCl and KCl-KI," *J. Phys. Chem.* **59**(2), 118-127 (1955). [DOI](https://doi.org/10.1021/j150524a007)

[3] P. Fellner, I. Votava, and M. Chrenkova-Paucirova, "Density of the Molten System LiCl-NaCl-KCl," *Chemicke Zvesti / Chemical Papers* **34**(3), 330-334 (1980). [Article record](https://www.chempap.org/index.php/chempap/article/view/1898), [full paper](https://www.chempap.org/file_access.php?file=343a330.pdf).

[4] F. Kohler, "Zur Berechnung der thermodynamischen Daten eines ternaren Systems aus den zugehorigen binaren Systemen," *Monatsh. Chem.* **91**, 738-740 (1960). [DOI](https://doi.org/10.1007/BF00899814)

[5] O. Redlich and A. T. Kister, "Algebraic Representation of Thermodynamic Properties and the Classification of Solutions," *Ind. Eng. Chem.* **40**(2), 345-348 (1948). [DOI](https://doi.org/10.1021/ie50458a036)

[6] M. J. Mehl, D. Hicks, C. Toher, O. Levy, R. M. Hanson, G. Hart, and S. Curtarolo, "The AFLOW Library of Crystallographic Prototypes: Part 1," *Comput. Mater. Sci.* **136**, S1-S828 (2017). [DOI](https://doi.org/10.1016/j.commatsci.2017.01.017), [rock-salt prototype](https://aflow.org/p/AB_cF8_225_a_b-001/).

[7] VASP Wiki, "POSCAR." [Format documentation](https://vasp.at/wiki/POSCAR)

[8] S. R. Hall, F. H. Allen, and I. D. Brown, "The Crystallographic Information File (CIF): a New Standard Archive File for Crystallography," *Acta Cryst. A* **47**, 655-685 (1991). [DOI](https://doi.org/10.1107/S010876739101067X)

[9] LAMMPS documentation: [`atom_style`](https://docs.lammps.org/atom_style.html), [`read_data`](https://docs.lammps.org/read_data.html).

[10] Bureau International des Poids et Mesures, *The International System of Units (SI Brochure)*, 9th ed. [DOI](https://doi.org/10.59161/AUEZ1291)

[11] NIST Chemistry WebBook, SRD 69. [DOI](https://doi.org/10.18434/T4D303), [WebBook](https://webbook.nist.gov/chemistry/).

[12] G. J. Janz, "Thermodynamic and Transport Properties for Molten Salts: Correlation Equations for Critically Evaluated Density, Surface Tension, Electrical Conductance, and Viscosity Data," *J. Phys. Chem. Ref. Data* **17**, Supplement 2 (1988). [NIST supplement index](https://www.nist.gov/srd/journal-physical-and-chemical-reference-data-monographs-or-supplements), [PDF](https://srd.nist.gov/JPCRD/jpcrdS2Vol17.pdf).

\normalsize
