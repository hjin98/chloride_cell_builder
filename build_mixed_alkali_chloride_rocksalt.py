#!/usr/bin/env python3
"""Generate mixed group-I/group-II chloride initial melt cells.

Supported cations are the stable alkali/alkaline-earth metals for which the
packaged molten-chloride volume model has sourced density data:
Li+, Na+, K+, Rb+, Cs+, Be2+, Mg2+, Ca2+, Sr2+, and Ba2+.

The cation lattice contains N = 4*nx*ny*nz sites.  Chloride count is not
hard-coded to N: charge neutrality is enforced exactly as

    N_Cl = sum_i z_i N_i,

where z_i is +1 or +2.  Monovalent-only cells reproduce the original rocksalt
site topology.  When divalent cations are present, the additional chloride ions
are placed on a second, tetrahedral-like auxiliary sublattice.  The result is a
compact crystalline starting template for subsequent melting/equilibration,
not a claim about the equilibrium solid structure of the mixture.

Volume model
------------
For every supported chloride component i (ACl or BCl2), the pure molten molar
volume is

    V_i^0(T) = M_i / rho_i^0(T).

The general composition model interpolates these pure molar volumes on a
per-cation (equivalently, per salt formula-unit) mole-fraction basis:

    V_m^ideal(T,x) = sum_i x_i V_i^0(T).

The previously implemented LiCl-NaCl-KCl excess-volume treatment is retained
exactly as a legacy compatibility model: source binary density correlations are
converted to binary excess molar volume and combined with the historical
symmetric Kohler-style geometric weighting.  The assembled ternary interior is
a package interpolation; it is not the direct Fellner et al. ternary equation.
This distinction is deliberate because prior calculations may depend on the
legacy numerical field.  In mixtures that also contain other cations, those
same known Li/Na/K pair contributions are retained with their existing
pair-sum-squared weights.  No unsupported excess-volume interaction is invented
for Rb/Cs/alkaline-earth pairs.

Density sources
---------------
[R1,R2,R4,R5] document the legacy LiCl-NaCl-KCl input equations and
interpolation-method provenance. [R3] is an external ternary validation/comparison
source; the default ternary path does not implement the Fellner equation.
[R12] supplies critically evaluated pure-melt density correlations for the
new chlorides.  The compact tabulation used here expresses them as
rho(t_C) = rho_m - k*(t_C-t_m), with the published validity ranges.

The radioactive Fr+/Ra2+ endpoints and H+ are intentionally not assigned an
empirical molten-chloride volume model: no compatible critically evaluated
pure-melt density correlation is provided by the source set used here.
"""

# SSDP 6.6.0 authority map for this executable:
#   D1 scientific meaning        -> SCIENTIFIC_METHOD.md
#   D2 numerical method          -> NUMERICAL_METHOD.md
#   D3 architecture              -> ARCHITECTURE.md
#   D4 public software contract  -> SPECIFICATION.md
#   external source locators     -> REFERENCES.md
# The Python implementation concretizes those owners; code/tests are evidence,
# not a substitute for the upstream scientific or numerical authority.

from __future__ import annotations

import argparse
import math
import random
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

__version__ = "1.3.4"

AVOGADRO = 6.02214076e23
CM3_TO_A3 = 1.0e24

# Keep the legacy Li/Na/K/Cl masses unchanged so old cells remain numerically
# reproducible.  New entries use standard atomic weights compatible with the
# same precision level.
ATOMIC_MASS = {
    "Li": 6.94100,
    "Na": 22.98977,
    "K": 39.09830,
    "Rb": 85.46780,
    "Cs": 132.90545196,
    "Be": 9.0121831,
    "Mg": 24.30500,
    "Ca": 40.07800,
    "Sr": 87.62000,
    "Ba": 137.32700,
    "Cl": 35.45300,
}

CATION_CHARGE = {
    "Li": 1,
    "Na": 1,
    "K": 1,
    "Rb": 1,
    "Cs": 1,
    "Be": 2,
    "Mg": 2,
    "Ca": 2,
    "Sr": 2,
    "Ba": 2,
}
CATION_ORDER = tuple(CATION_CHARGE)
CHARGES = {**{s: float(z) for s, z in CATION_CHARGE.items()}, "Cl": -1.0}

SALT_FORMULA = {s: f"{s}Cl" if z == 1 else f"{s}Cl2" for s, z in CATION_CHARGE.items()}
MOLAR_MASS = {
    SALT_FORMULA[s]: ATOMIC_MASS[s] + CATION_CHARGE[s] * ATOMIC_MASS["Cl"]
    for s in CATION_ORDER
}
# Exact legacy formula masses from v1.2.0.
MOLAR_MASS.update({"LiCl": 42.39400, "NaCl": 58.44277, "KCl": 74.55130})

# Source-support ranges for the legacy Li/Na/K correlations [R1].
#
# Pure LiCl and NaCl are the actual Table 500 endpoint tabulations used by the
# corresponding linear correlations.  Pure KCl is evaluated from the NaCl-KCl
# two-variable polynomial at C=100 mol% KCl, so its endpoint range is taken
# from that source table rather than from a different KCl correlation.
LEGACY_PURE_T_RANGES_K = {
    "LiCl": (880.0, 1070.0),
    "NaCl": (1090.0, 1220.0),
    "KCl": (1060.0, 1200.0),
}

# The LiCl-NaCl arbitrary-composition model is constructed from the 50:50
# source correlation, whose reported study range is 893-1123 K.  The NaCl-KCl
# and KCl-LiCl entries are broad source-study envelopes for the fitted binary
# surfaces; their exact liquid-data interval is composition dependent.
PAIR_T_RANGES_K = {
    ("LiCl", "NaCl"): (893.0, 1123.0),
    ("NaCl", "KCl"): (943.0, 1300.0),
    ("KCl", "LiCl"): (668.0, 1212.0),
}

# Prevent binary floating-point conversion (e.g. K -> C at an exact tabulated
# endpoint) from producing a false extrapolation warning.  This is far below
# any physically meaningful temperature precision.
TEMPERATURE_RANGE_ATOL = 1.0e-9

# Janz 1988 pure-melt density table [R12], expressed as
# rho(T_C) = rho_m - k*(T_C - T_m_C).  The model uses these entries for the
# newly added species.  Li/Na/K keep their legacy v1.2.0 correlations.
# Values are (T_m_C, rho_m_g_cm3, k_g_cm3_C, T_max_C).
PURE_DENSITY_LINEAR_C = {
    "RbCl": (715.0, 2.248, 0.000883, 923.0),
    "CsCl": (645.0, 2.790, 0.001065, 906.0),
    "BeCl2": (415.0, 1.540, 0.001100, 473.0),
    "MgCl2": (714.0, 1.680, 0.000271, 826.0),
    "CaCl2": (775.0, 2.085, 0.000422, 950.0),
    "SrCl2": (874.0, 2.727, 0.000578, 1037.0),
    "BaCl2": (962.0, 3.174, 0.000681, 1081.0),
}

CL_FRAC_CONV = [
    (0.0, 0.0, 0.0),
    (0.0, 0.5, 0.5),
    (0.5, 0.0, 0.5),
    (0.5, 0.5, 0.0),
]
M_FRAC_CONV = [
    (0.5, 0.0, 0.0),
    (0.0, 0.5, 0.0),
    (0.0, 0.0, 0.5),
    (0.5, 0.5, 0.5),
]
# One auxiliary chloride site per cation site.  Occupying all of them gives the
# extra N chlorides needed for a purely divalent composition.
EXTRA_CL_FRAC_CONV = [
    (0.25, 0.25, 0.25),
    (0.25, 0.75, 0.75),
    (0.75, 0.25, 0.75),
    (0.75, 0.75, 0.25),
]


@dataclass
class Atom:
    symbol: str
    frac: Tuple[float, float, float]


# -------------------------- density equations ---------------------------

def rho_pure_LiCl(T_K: float) -> float:
    return 1.8721 - 4.243e-4 * T_K


def rho_pure_NaCl(T_K: float) -> float:
    return 2.1390 - 5.444e-4 * T_K


def rho_pure_KCl(T_K: float) -> float:
    C = 100.0
    return (
        2.14250
        - 5.45436e-4 * T_K
        + 2.79726e-8 * C**3
        - 9.89901e-7 * T_K * C
        + 2.78444e-10 * C * T_K**2
    )


def rho_pure_extended(salt: str, T_K: float) -> float:
    """Pure molten density from the Janz 1988 linear-in-Celsius tabulation [R12]."""
    try:
        t_m_c, rho_m, k, _t_max_c = PURE_DENSITY_LINEAR_C[salt]
    except KeyError as exc:
        raise ValueError(f"No extended pure density correlation for {salt}") from exc
    T_C = T_K - 273.15
    return rho_m - k * (T_C - t_m_c)


def rho_pure_salt(T_K: float, salt: str) -> float:
    if salt == "LiCl":
        return rho_pure_LiCl(T_K)
    if salt == "NaCl":
        return rho_pure_NaCl(T_K)
    if salt == "KCl":
        return rho_pure_KCl(T_K)
    return rho_pure_extended(salt, T_K)


def rho_binary_LiCl_NaCl(T_K: float, y_na: float) -> float:
    y_li = 1.0 - y_na
    V_ideal = y_li * MOLAR_MASS["LiCl"] / rho_pure_LiCl(T_K) + y_na * MOLAR_MASS["NaCl"] / rho_pure_NaCl(T_K)
    rho_mid = 2.0121 - 4.858e-4 * T_K
    V_mid = (0.5 * MOLAR_MASS["LiCl"] + 0.5 * MOLAR_MASS["NaCl"]) / rho_mid
    V_mid_ideal = 0.5 * MOLAR_MASS["LiCl"] / rho_pure_LiCl(T_K) + 0.5 * MOLAR_MASS["NaCl"] / rho_pure_NaCl(T_K)
    L0 = 4.0 * (V_mid - V_mid_ideal)
    V = V_ideal + y_li * y_na * L0
    return (y_li * MOLAR_MASS["LiCl"] + y_na * MOLAR_MASS["NaCl"]) / V


def rho_binary_NaCl_KCl(T_K: float, y_k: float) -> float:
    C = 100.0 * y_k
    return (
        2.14250
        - 5.45436e-4 * T_K
        + 2.79726e-8 * C**3
        - 9.89901e-7 * T_K * C
        + 2.78444e-10 * C * T_K**2
    )


def rho_binary_KCl_LiCl(T_K: float, y_li: float) -> float:
    C = 100.0 * y_li
    return (
        1.98715
        - 5.89069e-4 * T_K
        + 2.02443e-3 * C
        + 5.31098e-11 * T_K**3
        + 1.21294e-7 * C**3
        - 2.16826e-8 * T_K * C**2
    )


def pure_molar_volume_cm3_mol(salt: str, T_K: float) -> float:
    rho = rho_pure_salt(T_K, salt)
    if rho <= 0:
        raise ValueError(f"Non-positive pure density for {salt} at T={T_K:g} K: rho={rho:g}")
    return MOLAR_MASS[salt] / rho


def binary_excess_from_density(T_K: float, salt_i: str, salt_j: str, y_i: float, y_j: float) -> float:
    if y_i <= 0.0 or y_j <= 0.0:
        return 0.0
    if (salt_i, salt_j) == ("LiCl", "NaCl"):
        rho = rho_binary_LiCl_NaCl(T_K, y_j)
    elif (salt_i, salt_j) == ("NaCl", "KCl"):
        rho = rho_binary_NaCl_KCl(T_K, y_j)
    elif (salt_i, salt_j) == ("KCl", "LiCl"):
        rho = rho_binary_KCl_LiCl(T_K, y_j)
    else:
        raise ValueError(f"Unsupported binary pair orientation: {salt_i}-{salt_j}")
    V_bin = (y_i * MOLAR_MASS[salt_i] + y_j * MOLAR_MASS[salt_j]) / rho
    V_ideal = y_i * pure_molar_volume_cm3_mol(salt_i, T_K) + y_j * pure_molar_volume_cm3_mol(salt_j, T_K)
    return V_bin - V_ideal


def mixture_molar_volume_from_fractions(T_K: float, cation_fractions: Mapping[str, float], use_excess: bool = True) -> float:
    """Mixture molar volume per mole of cations/salt formula units.

    Composition interpolation is performed over pure ACl/BCl2 molar volumes.
    The existing LiCl-NaCl-KCl pair excess contributions are then added with
    the frozen legacy pair-sum-squared/Kohler-style construction.  This is a
    compatibility model, not a claim of identity with the Fellner ternary
    equation.  Other pairs remain ideal because no source-backed excess model
    is packaged for them.
    """
    x_salt = {SALT_FORMULA[s]: float(cation_fractions.get(s, 0.0)) for s in CATION_ORDER}
    V = sum(x * pure_molar_volume_cm3_mol(salt, T_K) for salt, x in x_salt.items() if x > 0.0)
    if not use_excess:
        return V
    for salt_i, salt_j in [("LiCl", "NaCl"), ("NaCl", "KCl"), ("KCl", "LiCl")]:
        pair_sum = x_salt[salt_i] + x_salt[salt_j]
        if pair_sum <= 1e-15:
            continue
        y_i = x_salt[salt_i] / pair_sum
        y_j = x_salt[salt_j] / pair_sum
        V += pair_sum**2 * binary_excess_from_density(T_K, salt_i, salt_j, y_i, y_j)
    return V


def mixture_molar_volume_cm3_mol(T_K: float, x_li: float, x_na: float, x_k: float, use_excess: bool = True) -> float:
    """Backward-compatible Li/Na/K wrapper retained from v1.2.0."""
    return mixture_molar_volume_from_fractions(
        T_K,
        {"Li": x_li, "Na": x_na, "K": x_k},
        use_excess=use_excess,
    )


def direct_binary_density_if_edge(T_K: float, x_li: float, x_na: float, x_k: float) -> Optional[float]:
    eps = 1e-12
    if x_k < eps and x_li > eps and x_na > eps:
        return rho_binary_LiCl_NaCl(T_K, x_na / (x_li + x_na))
    if x_li < eps and x_na > eps and x_k > eps:
        return rho_binary_NaCl_KCl(T_K, x_k / (x_na + x_k))
    if x_na < eps and x_k > eps and x_li > eps:
        return rho_binary_KCl_LiCl(T_K, x_li / (x_k + x_li))
    return None


def direct_known_binary_density(T_K: float, cation_fractions: Mapping[str, float]) -> Optional[Tuple[str, float]]:
    active = [s for s in CATION_ORDER if cation_fractions.get(s, 0.0) > 1e-12]
    if len(active) != 2 or not set(active).issubset({"Li", "Na", "K"}):
        return None
    x_li = cation_fractions.get("Li", 0.0)
    x_na = cation_fractions.get("Na", 0.0)
    x_k = cation_fractions.get("K", 0.0)
    rho = direct_binary_density_if_edge(T_K, x_li, x_na, x_k)
    if rho is None:
        return None
    return "-".join(SALT_FORMULA[s] for s in active), rho


def mixture_molar_mass_from_fractions(cation_fractions: Mapping[str, float]) -> float:
    return sum(
        float(cation_fractions.get(s, 0.0)) * MOLAR_MASS[SALT_FORMULA[s]]
        for s in CATION_ORDER
    )


def mixture_molar_mass_g_mol(x_li: float, x_na: float, x_k: float) -> float:
    """Backward-compatible Li/Na/K wrapper retained from v1.2.0."""
    return mixture_molar_mass_from_fractions({"Li": x_li, "Na": x_na, "K": x_k})


# -------------------------- composition parsing --------------------------

def _provided_count(vals: Iterable[Optional[float]]) -> int:
    return sum(v is not None for v in vals)


def complete_fractions(x_li: Optional[float], x_na: Optional[float], x_k: Optional[float]) -> Tuple[float, float, float]:
    vals = [x_li, x_na, x_k]
    n_given = _provided_count(vals)
    if n_given == 1:
        only = next(v for v in vals if v is not None)
        if abs(float(only) - 1.0) > 1e-8:
            raise ValueError("A single supplied Li/Na/K mole fraction must equal 1.")
        vals = [0.0 if v is None else v for v in vals]
    elif n_given < 1:
        raise ValueError("Specify at least one mole fraction, or use explicit atom counts.")
    elif n_given == 2:
        missing = vals.index(None)
        vals[missing] = 1.0 - sum(v for v in vals if v is not None)
    vals_f = [float(v) for v in vals]  # type: ignore[arg-type]
    if any(v < -1e-12 for v in vals_f):
        raise ValueError(f"Mole fractions must be non-negative after completion: {vals_f}")
    total = sum(vals_f)
    if abs(total - 1.0) > 1e-8:
        raise ValueError(f"Mole fractions must sum to 1. Got {total:.16g} from {vals_f}.")
    vals_f = [0.0 if abs(v) < 1e-14 else v for v in vals_f]
    return vals_f[0], vals_f[1], vals_f[2]


def normalize_fraction_map(raw: Mapping[str, Optional[float]]) -> Dict[str, float]:
    supplied = {s: float(v) for s, v in raw.items() if v is not None}
    if not supplied:
        raise ValueError("Specify at least one cation mole fraction.")
    if any(not math.isfinite(v) or v < 0.0 for v in supplied.values()):
        raise ValueError(f"Mole fractions must be finite and non-negative: {supplied}")

    has_extended = any(s not in {"Li", "Na", "K"} for s in supplied)
    if not has_extended and len(supplied) <= 3:
        x_li, x_na, x_k = complete_fractions(raw.get("Li"), raw.get("Na"), raw.get("K"))
        result = {s: 0.0 for s in CATION_ORDER}
        result.update({"Li": x_li, "Na": x_na, "K": x_k})
        return result

    total = sum(supplied.values())
    if abs(total - 1.0) > 1e-8:
        raise ValueError(
            "When any extended cation is used, omitted species are zero and the supplied cation mole fractions "
            f"must sum to 1. Got {total:.16g}."
        )
    return {s: (0.0 if abs(supplied.get(s, 0.0)) < 1e-14 else supplied.get(s, 0.0)) for s in CATION_ORDER}


def counts_from_fraction_map(N: int, fractions: Mapping[str, float]) -> Dict[str, int]:
    raw = [N * float(fractions.get(s, 0.0)) for s in CATION_ORDER]
    base = [int(math.floor(v)) for v in raw]
    remainder = N - sum(base)
    order = sorted(range(len(CATION_ORDER)), key=lambda i: (raw[i] - base[i], -i), reverse=True)
    for i in order[:remainder]:
        base[i] += 1
    return dict(zip(CATION_ORDER, base))


def counts_from_fractions(N: int, x_li: float, x_na: float, x_k: float) -> Tuple[int, int, int]:
    counts = counts_from_fraction_map(N, {"Li": x_li, "Na": x_na, "K": x_k})
    return counts["Li"], counts["Na"], counts["K"]


def complete_counts(N: int, n_li: Optional[int], n_na: Optional[int], n_k: Optional[int]) -> Tuple[int, int, int]:
    vals = [n_li, n_na, n_k]
    n_given = _provided_count(vals)
    if n_given == 1:
        only = next(v for v in vals if v is not None)
        if int(only) != N:
            raise ValueError(f"A single supplied Li/Na/K count must equal N={N}.")
        vals = [0 if v is None else v for v in vals]
    elif n_given < 1:
        raise ValueError("Specify at least one explicit count, or use mole fractions.")
    elif n_given == 2:
        missing = vals.index(None)
        vals[missing] = N - sum(v for v in vals if v is not None)
    vals_i = [int(v) for v in vals]  # type: ignore[arg-type]
    if any(v < 0 for v in vals_i):
        raise ValueError(f"Atom counts must be non-negative after completion: {vals_i}")
    if sum(vals_i) != N:
        raise ValueError(f"Cation counts must sum to N={N}. Got {sum(vals_i)} from {vals_i}.")
    return vals_i[0], vals_i[1], vals_i[2]


def normalize_count_map(N: int, raw: Mapping[str, Optional[int]]) -> Dict[str, int]:
    supplied = {s: int(v) for s, v in raw.items() if v is not None}
    if not supplied:
        raise ValueError("Specify at least one explicit cation count.")
    if any(v < 0 for v in supplied.values()):
        raise ValueError(f"Cation counts must be non-negative: {supplied}")

    has_extended = any(s not in {"Li", "Na", "K"} for s in supplied)
    if not has_extended and len(supplied) <= 3:
        n_li, n_na, n_k = complete_counts(N, raw.get("Li"), raw.get("Na"), raw.get("K"))
        result = {s: 0 for s in CATION_ORDER}
        result.update({"Li": n_li, "Na": n_na, "K": n_k})
        return result

    total = sum(supplied.values())
    if total != N:
        raise ValueError(
            "When any extended cation is used, omitted species are zero and the supplied cation counts "
            f"must sum to N={N}. Got {total}."
        )
    return {s: supplied.get(s, 0) for s in CATION_ORDER}


def fractions_from_counts(counts: Mapping[str, int], N: int) -> Dict[str, float]:
    return {s: counts.get(s, 0) / N for s in CATION_ORDER}


def required_chloride_count(counts: Mapping[str, int]) -> int:
    return sum(CATION_CHARGE[s] * int(counts.get(s, 0)) for s in CATION_ORDER)


# ---------------------------- structure build ----------------------------

def build_rocksalt_atoms(nx: int, ny: int, nz: int, cation_symbols: Sequence[str]) -> List[Atom]:
    expected = 4 * nx * ny * nz
    if len(cation_symbols) != expected:
        raise ValueError(f"Expected {expected} cations, got {len(cation_symbols)}")
    atoms: List[Atom] = []
    cation_index = 0
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                for fx, fy, fz in CL_FRAC_CONV:
                    atoms.append(Atom("Cl", ((ix + fx) / nx, (iy + fy) / ny, (iz + fz) / nz)))
                for fx, fy, fz in M_FRAC_CONV:
                    atoms.append(Atom(cation_symbols[cation_index], ((ix + fx) / nx, (iy + fy) / ny, (iz + fz) / nz)))
                    cation_index += 1
    return atoms


def build_mixed_chloride_atoms(
    nx: int,
    ny: int,
    nz: int,
    cation_symbols: Sequence[str],
    rng: random.Random,
) -> List[Atom]:
    """Build the legacy rocksalt template plus charge-required extra Cl sites."""
    atoms = build_rocksalt_atoms(nx, ny, nz, cation_symbols)
    extra_count = sum(CATION_CHARGE[s] - 1 for s in cation_symbols)
    if extra_count == 0:
        return atoms

    candidates: List[Atom] = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                for fx, fy, fz in EXTRA_CL_FRAC_CONV:
                    candidates.append(Atom("Cl", ((ix + fx) / nx, (iy + fy) / ny, (iz + fz) / nz)))
    if extra_count > len(candidates):
        raise ValueError(f"Need {extra_count} auxiliary chloride sites but only {len(candidates)} exist.")
    for idx in rng.sample(range(len(candidates)), extra_count):
        atoms.append(candidates[idx])
    return atoms


def lattice_vectors(a_conv_A: float, nx: int, ny: int, nz: int) -> Tuple[Tuple[float, float, float], ...]:
    return ((a_conv_A * nx, 0.0, 0.0), (0.0, a_conv_A * ny, 0.0), (0.0, 0.0, a_conv_A * nz))


def frac_to_cart(frac: Tuple[float, float, float], lattice: Tuple[Tuple[float, float, float], ...]) -> Tuple[float, float, float]:
    return (frac[0] * lattice[0][0], frac[1] * lattice[1][1], frac[2] * lattice[2][2])


# ------------------------------- writers ---------------------------------

def species_order(atoms: Sequence[Atom]) -> List[str]:
    preferred = list(CATION_ORDER) + ["Cl"]
    present = {a.symbol for a in atoms}
    return [s for s in preferred if s in present] + sorted(present - set(preferred))


def write_poscar(path: Path, atoms: Sequence[Atom], lattice: Tuple[Tuple[float, float, float], ...], comment: str) -> None:
    order = species_order(atoms)
    grouped = [a for s in order for a in atoms if a.symbol == s]
    with path.open("w") as f:
        f.write(comment[:240] + "\n")
        f.write("1.0\n")
        for v in lattice:
            f.write(f"{v[0]:20.12f} {v[1]:20.12f} {v[2]:20.12f}\n")
        f.write(" ".join(order) + "\n")
        f.write(" ".join(str(sum(a.symbol == s for a in atoms)) for s in order) + "\n")
        f.write("Direct\n")
        for a in grouped:
            f.write(f"{a.frac[0]:18.12f} {a.frac[1]:18.12f} {a.frac[2]:18.12f}\n")


def write_cif(path: Path, atoms: Sequence[Atom], lattice: Tuple[Tuple[float, float, float], ...], comment: str) -> None:
    a, b, c = lattice[0][0], lattice[1][1], lattice[2][2]
    with path.open("w") as f:
        f.write("data_mixed_group12_chloride_template\n")
        f.write(f"_audit_creation_method '{comment}'\n")
        f.write("_symmetry_space_group_name_H-M 'P 1'\n")
        f.write("_symmetry_Int_Tables_number 1\n")
        f.write(f"_cell_length_a {a:.12f}\n")
        f.write(f"_cell_length_b {b:.12f}\n")
        f.write(f"_cell_length_c {c:.12f}\n")
        f.write("_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n")
        f.write("loop_\n_symmetry_equiv_pos_as_xyz\n'x, y, z'\n")
        f.write("loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_occupancy\n")
        counters: Dict[str, int] = {}
        for atom in atoms:
            counters[atom.symbol] = counters.get(atom.symbol, 0) + 1
            label = f"{atom.symbol}{counters[atom.symbol]}"
            f.write(f"{label:8s} {atom.symbol:2s} {atom.frac[0]:.12f} {atom.frac[1]:.12f} {atom.frac[2]:.12f} 1.0\n")


def write_lammps_data(
    path: Path,
    atoms: Sequence[Atom],
    lattice: Tuple[Tuple[float, float, float], ...],
    comment: str,
    atom_style: str,
) -> None:
    order = species_order(atoms)
    type_id = {s: i + 1 for i, s in enumerate(order)}
    lx, ly, lz = lattice[0][0], lattice[1][1], lattice[2][2]
    with path.open("w") as f:
        f.write(comment + "\n\n")
        f.write(f"{len(atoms)} atoms\n")
        f.write(f"{len(order)} atom types\n\n")
        f.write(f"0.0 {lx:.12f} xlo xhi\n")
        f.write(f"0.0 {ly:.12f} ylo yhi\n")
        f.write(f"0.0 {lz:.12f} zlo zhi\n\n")
        f.write("Masses\n\n")
        for s in order:
            f.write(f"{type_id[s]} {ATOMIC_MASS[s]:.8f} # {s}\n")
        f.write(f"\nAtoms # {atom_style}\n\n")
        for i, atom in enumerate(atoms, start=1):
            x, y, z = frac_to_cart(atom.frac, lattice)
            if atom_style == "charge":
                f.write(f"{i} {type_id[atom.symbol]} {CHARGES[atom.symbol]:.6f} {x:.12f} {y:.12f} {z:.12f} # {atom.symbol}\n")
            elif atom_style == "atomic":
                f.write(f"{i} {type_id[atom.symbol]} {x:.12f} {y:.12f} {z:.12f} # {atom.symbol}\n")
            else:
                raise ValueError(f"Unsupported LAMMPS atom style: {atom_style}")


# ----------------------------- CLI/helpers -------------------------------

def parse_hydrostatic_volume_strain(value: str) -> float:
    raw = value.strip()
    is_percent = raw.endswith("%")
    number_text = raw[:-1].strip() if is_percent else raw
    if not number_text:
        raise argparse.ArgumentTypeError("Hydrostatic strain must contain a numeric value.")
    try:
        strain = float(number_text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"Invalid hydrostatic strain {value!r}; use a fraction such as 0.05 or a percentage such as +5%."
        ) from exc
    if is_percent:
        strain /= 100.0
    if not math.isfinite(strain):
        raise argparse.ArgumentTypeError("Hydrostatic strain must be finite.")
    if strain <= -1.0:
        raise argparse.ArgumentTypeError(
            "Hydrostatic volume strain must be greater than -1 (-100%) so the strained volume remains positive."
        )
    return strain


def _normalize_strain_argv(argv: Sequence[str]) -> List[str]:
    flags = {"--hydrostatic-strain", "--volume-strain"}
    normalized: List[str] = []
    i = 0
    while i < len(argv):
        token = argv[i]
        if token in flags and i + 1 < len(argv):
            next_token = argv[i + 1].strip()
            if next_token.startswith("-") and next_token.endswith("%"):
                normalized.append(f"{token}={next_token}")
                i += 2
                continue
        normalized.append(token)
        i += 1
    return normalized


def apply_hydrostatic_volume_strain(volume: float, volume_strain: float) -> Tuple[float, float]:
    if not math.isfinite(volume_strain) or volume_strain <= -1.0:
        raise ValueError("Hydrostatic volume strain must be finite and greater than -1.")
    volume_factor = 1.0 + volume_strain
    linear_factor = volume_factor ** (1.0 / 3.0)
    return volume * volume_factor, linear_factor


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Build a randomized Li/Na/K/Rb/Cs/Be/Mg/Ca/Sr/Ba chloride initial melt template, "
            "charge balanced with Cl- and scaled to a temperature/composition-dependent molten volume."
        ),
    )
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    p.add_argument("--nx", type=int, required=True, help="Number of conventional template cells along a.")
    p.add_argument("--ny", type=int, required=True, help="Number of conventional template cells along b.")
    p.add_argument("--nz", type=int, required=True, help="Number of conventional template cells along c.")
    tgroup = p.add_mutually_exclusive_group(required=True)
    tgroup.add_argument("--temp-c", type=float, help="Target liquid temperature in degrees Celsius.")
    tgroup.add_argument("--temp-k", type=float, help="Target liquid temperature in Kelvin.")

    comp = p.add_argument_group(
        "composition by cation mole fraction; omitted extended species are zero. "
        "Legacy Li/Na/K two-value completion is retained"
    )
    for s in CATION_ORDER:
        comp.add_argument(f"--x-{s.lower()}", type=float, default=None, help=f"Cation mole fraction of {SALT_FORMULA[s]} component.")

    cnt = p.add_argument_group(
        "composition by explicit cation count; counts refer to cations, while Cl- is generated automatically for charge balance"
    )
    for s in CATION_ORDER:
        cnt.add_argument(f"--n-{s.lower()}", type=int, default=None, help=f"Number of {s} cations.")

    p.add_argument("--seed", type=int, default=12345, help="Random seed for cation and auxiliary-Cl site population.")
    p.add_argument("--format", choices=["poscar", "cif", "lammps", "all"], default="all", help="Output structural format.")
    p.add_argument("--output-prefix", default="mixed_group12_chloride", help="Output filename prefix. A final temperature tag is appended/replaced unless --no-auto-temp-tag is used.")
    p.add_argument("--no-auto-temp-tag", action="store_true", help="Do not append or replace the temperature tag in the output prefix.")
    p.add_argument(
        "--ideal-only",
        action="store_true",
        help="Use pure-component molar-volume interpolation only; omit the source-backed Li/Na/K excess-volume terms.",
    )
    p.add_argument(
        "--hydrostatic-strain",
        "--volume-strain",
        dest="hydrostatic_strain",
        type=parse_hydrostatic_volume_strain,
        default=0.0,
        metavar="EPS_V",
        help=(
            "Hydrostatic strain applied after the equilibrium liquid volume is predicted, using a volume-strain "
            "convention. Bare values are fractional volume strain (0.05 = +5%% volume); percent-suffixed values "
            "are also accepted (+5%% = 0.05). Positive expands, negative compresses; value must be > -100%%."
        ),
    )
    p.add_argument("--lammps-atom-style", choices=["atomic", "charge"], default="atomic", help="LAMMPS Atoms-section format.")
    if argv is None:
        argv = sys.argv[1:]
    return p.parse_args(_normalize_strain_argv(list(argv)))


def sanitize_float_for_name(value: float) -> str:
    if abs(value - round(value)) < 1e-10:
        return str(int(round(value)))
    s = f"{value:.8g}"
    return s.replace("-", "m").replace(".", "p")


def temperature_values(args: argparse.Namespace) -> Tuple[float, float, str]:
    if args.temp_k is not None:
        T_K = float(args.temp_k)
        T_C = T_K - 273.15
        tag = f"{sanitize_float_for_name(T_K)}K"
    else:
        T_C = float(args.temp_c)
        T_K = T_C + 273.15
        tag = f"{sanitize_float_for_name(T_C)}C"
    if T_K <= 0:
        raise SystemExit(f"Temperature must be positive in Kelvin. Got {T_K:g} K.")
    return T_K, T_C, tag


def output_prefix_with_temperature(prefix: str, temp_tag: str, no_auto: bool) -> Path:
    if no_auto:
        return Path(prefix)
    new_prefix = re.sub(r"_(?:m?\d+(?:p\d+)?)(?:C|K)$", "", prefix)
    return Path(f"{new_prefix}_{temp_tag}")


def format_counts(counts: Mapping[str, int]) -> str:
    return ",".join(f"{s}:{counts[s]}" for s in CATION_ORDER if counts.get(s, 0) > 0)


def format_fractions(fractions: Mapping[str, float]) -> str:
    return ",".join(f"{s}:{fractions[s]:.12f}" for s in CATION_ORDER if fractions.get(s, 0.0) > 0.0)


def make_comment(
    T_K: float,
    T_C: float,
    nx: int,
    ny: int,
    nz: int,
    counts: Mapping[str, int],
    n_cl: int,
    Vm_eq: float,
    volume_strain: float,
    Vm_strained: float,
    lammps_atom_style: str,
) -> str:
    return (
        f"Group-I/II chloride initial melt template; nx,ny,nz={nx},{ny},{nz}; "
        f"T={T_K:g} K ({T_C:g} C); cations={format_counts(counts)}; Cl={n_cl}; "
        f"Vm_eq={Vm_eq:.8f} cm3/mol-cation; epsV={volume_strain:+.8g}; "
        f"Vm={Vm_strained:.8f} cm3/mol-cation; lammps_atom_style={lammps_atom_style}"
    )


def _outside_temperature_range(value: float, lo: float, hi: float) -> bool:
    """Return True only for a material excursion beyond an inclusive range."""
    return value < lo - TEMPERATURE_RANGE_ATOL or value > hi + TEMPERATURE_RANGE_ATOL


def warn_temperature_ranges(
    T_K: float,
    cation_fractions: Mapping[str, float],
    use_excess: bool = True,
) -> List[str]:
    """Report empirical-range excursions without changing the requested state point.

    The legacy Li/Na/K numerical field is deliberately frozen.  These warnings
    therefore diagnose extrapolation; they do not clamp temperature or alter
    the volume model.  Every active Li/Na/K pure-component correlation is
    checked because it contributes to the ideal-volume baseline, and every
    active legacy binary excess correlation is checked separately when the
    excess model is enabled.
    """
    warnings: List[str] = []
    active = {s: float(cation_fractions.get(s, 0.0)) for s in CATION_ORDER}

    # Every active legacy pure correlation contributes to the ideal-volume
    # baseline, including inside mixtures, so diagnose its source range wherever
    # it is actually evaluated.
    for species in ("Li", "Na", "K"):
        if active[species] <= 1e-14:
            continue
        salt = SALT_FORMULA[species]
        lo, hi = LEGACY_PURE_T_RANGES_K[salt]
        if _outside_temperature_range(T_K, lo, hi):
            warnings.append(
                f"T={T_K:g} K is outside the Janz/NIST pure-{salt} density range "
                f"[{lo:g}, {hi:g}] K; using the legacy pure-component correlation as an extrapolation."
            )

    if use_excess:
        x_li = active["Li"]
        x_na = active["Na"]
        x_k = active["K"]
        pair_data = [
            (("LiCl", "NaCl"), x_li * x_na),
            (("NaCl", "KCl"), x_na * x_k),
            (("KCl", "LiCl"), x_k * x_li),
        ]
        for pair, weight in pair_data:
            if weight <= 1e-14:
                continue
            lo, hi = PAIR_T_RANGES_K[pair]
            if _outside_temperature_range(T_K, lo, hi):
                if pair == ("LiCl", "NaCl"):
                    detail = "50:50 source-correlation range used to construct the legacy excess term"
                else:
                    detail = "broad source-study envelope; the exact liquid-data interval is composition-dependent"
                warnings.append(
                    f"T={T_K:g} K is outside the Janz/NIST {pair[0]}-{pair[1]} "
                    f"{detail} [{lo:g}, {hi:g}] K; using the legacy binary model as an extrapolation."
                )

    T_C = T_K - 273.15
    for s in CATION_ORDER:
        if active[s] <= 1e-14:
            continue
        salt = SALT_FORMULA[s]
        if salt in PURE_DENSITY_LINEAR_C:
            t_m_c, _rho_m, _k, t_max_c = PURE_DENSITY_LINEAR_C[salt]
            if _outside_temperature_range(T_C, t_m_c, t_max_c):
                warnings.append(
                    f"T={T_C:g} C is outside the Janz 1988 pure-{salt} density range "
                    f"[{t_m_c:g}, {t_max_c:g}] C; using the tabulated linear correlation as an extrapolation."
                )
    return warnings


def _arg_map(args: argparse.Namespace, prefix: str) -> Dict[str, Optional[float]]:
    return {s: getattr(args, f"{prefix}_{s.lower()}") for s in CATION_ORDER}


def main() -> None:
    args = parse_args()
    if args.nx <= 0 or args.ny <= 0 or args.nz <= 0:
        raise SystemExit("nx, ny, and nz must be positive integers.")

    T_K, T_C, temp_tag = temperature_values(args)
    N = 4 * args.nx * args.ny * args.nz

    fraction_inputs = _arg_map(args, "x")
    count_inputs = _arg_map(args, "n")
    has_fraction = any(v is not None for v in fraction_inputs.values())
    has_count = any(v is not None for v in count_inputs.values())
    if has_fraction and has_count:
        raise SystemExit("Use either cation mole fractions or explicit cation counts, not both.")
    if not has_fraction and not has_count:
        raise SystemExit("Specify composition using either cation mole fractions or explicit cation counts.")

    try:
        if has_count:
            counts = normalize_count_map(N, count_inputs)  # type: ignore[arg-type]
            fractions = fractions_from_counts(counts, N)
            input_mode = "explicit cation counts"
        else:
            requested_fractions = normalize_fraction_map(fraction_inputs)
            counts = counts_from_fraction_map(N, requested_fractions)
            fractions = fractions_from_counts(counts, N)
            input_mode = "cation mole fractions rounded to integer cation counts"
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc

    n_cl = required_chloride_count(counts)
    positive_charge = sum(CATION_CHARGE[s] * counts[s] for s in CATION_ORDER)
    if positive_charge != n_cl:
        raise SystemExit("Internal error: charge-neutral chloride count construction failed.")

    try:
        Vm_eq = mixture_molar_volume_from_fractions(T_K, fractions, use_excess=not args.ideal_only)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    Mmix = mixture_molar_mass_from_fractions(fractions)
    rho_eq = Mmix / Vm_eq
    equilibrium_volume_A3 = N * Vm_eq * CM3_TO_A3 / AVOGADRO

    Vm, linear_scale = apply_hydrostatic_volume_strain(Vm_eq, args.hydrostatic_strain)
    volume_factor = 1.0 + args.hydrostatic_strain
    rho = Mmix / Vm
    total_volume_A3 = equilibrium_volume_A3 * volume_factor
    n_conv = args.nx * args.ny * args.nz
    a_conv_A = (total_volume_A3 / n_conv) ** (1.0 / 3.0)
    lattice = lattice_vectors(a_conv_A, args.nx, args.ny, args.nz)

    cations = [s for s in CATION_ORDER for _ in range(counts[s])]
    rng = random.Random(args.seed)
    rng.shuffle(cations)
    atoms = build_mixed_chloride_atoms(args.nx, args.ny, args.nz, cations, rng)

    actual_cl = sum(a.symbol == "Cl" for a in atoms)
    if actual_cl != n_cl:
        raise SystemExit(f"Internal error: generated {actual_cl} Cl atoms but charge balance requires {n_cl}.")
    total_charge = sum(CHARGES[a.symbol] for a in atoms)
    if abs(total_charge) > 1e-12:
        raise SystemExit(f"Internal error: generated cell is not charge neutral (Q={total_charge:g} e).")

    prefix = output_prefix_with_temperature(args.output_prefix, temp_tag, args.no_auto_temp_tag)
    comment = make_comment(
        T_K,
        T_C,
        args.nx,
        args.ny,
        args.nz,
        counts,
        n_cl,
        Vm_eq,
        args.hydrostatic_strain,
        Vm,
        args.lammps_atom_style,
    )
    written: List[Path] = []
    if args.format in ("poscar", "all"):
        path = prefix.with_suffix(".POSCAR")
        write_poscar(path, atoms, lattice, comment)
        written.append(path)
    if args.format in ("cif", "all"):
        path = prefix.with_suffix(".cif")
        write_cif(path, atoms, lattice, comment)
        written.append(path)
    if args.format in ("lammps", "all"):
        path = prefix.with_suffix(".lammps.data")
        write_lammps_data(path, atoms, lattice, comment, args.lammps_atom_style)
        written.append(path)

    print("=== Mixed group-I/group-II chloride initial melt cell ===")
    print(f"Input mode                  : {input_mode}")
    print(f"Template cells nx ny nz     : {args.nx} {args.ny} {args.nz}")
    print(f"Cation sites N              : {N}")
    print(f"Chloride ions N_Cl          : {n_cl}")
    print(f"Atoms total                 : {N + n_cl}")
    print(f"Net ionic charge            : {total_charge:+.1f} e")
    print(f"Target temperature          : {T_K:.6f} K  ({T_C:.6f} C)")
    print(f"Cation counts               : {format_counts(counts)}")
    print(f"Achieved cation fractions   : {format_fractions(fractions)}")
    print(f"Equilibrium molar volume Vm0: {Vm_eq:.10f} cm^3/mol-cation")
    print(f"Mixture molar mass          : {Mmix:.10f} g/mol-cation-components")
    print(f"Equilibrium density rho0    : {rho_eq:.10f} g/cm^3")
    edge = direct_known_binary_density(T_K, fractions)
    if edge is not None and not args.ideal_only:
        pair_name, edge_rho = edge
        print(f"Direct active-binary rho0   : {edge_rho:.10f} g/cm^3 ({pair_name})")
    print(f"Hydrostatic volume strain   : {args.hydrostatic_strain:+.10f} ({100.0 * args.hydrostatic_strain:+.6f}%)")
    print(f"Volume scale factor         : {volume_factor:.10f}")
    print(f"Isotropic length factor     : {linear_scale:.10f}")
    print(f"Strained molar volume Vm    : {Vm:.10f} cm^3/mol-cation")
    print(f"Strained density rho        : {rho:.10f} g/cm^3")
    print(f"Equilibrium cell volume     : {equilibrium_volume_A3:.10f} A^3")
    print(f"Total cell volume           : {total_volume_A3:.10f} A^3")
    print(f"Equivalent template a       : {a_conv_A:.10f} A")
    print(f"Supercell lengths           : {lattice[0][0]:.10f} {lattice[1][1]:.10f} {lattice[2][2]:.10f} A")
    if args.ideal_only:
        print("Excess-volume terms         : off, pure-component molar-volume interpolation only")
    else:
        print("Excess-volume terms         : legacy Li/Na/K pairwise/Kohler compatibility model on; all other pairs ideal")
    if any(CATION_CHARGE[s] == 2 and counts[s] > 0 for s in CATION_ORDER):
        print("Template topology           : rocksalt base + auxiliary Cl sites required by +2 cations")
    else:
        print("Template topology           : legacy rocksalt")
    if args.format in ("lammps", "all"):
        print(f"LAMMPS atom style           : {args.lammps_atom_style}")
        print(f"LAMMPS input reminder       : use 'atom_style {args.lammps_atom_style}' before read_data")
    for warning in warn_temperature_ranges(T_K, fractions, use_excess=not args.ideal_only):
        print(f"WARNING: {warning}")
    print("Written files:")
    for path in written:
        print(f"  {path}")


# --------------------------- bibliographic references --------------------
# [R1] G. J. Janz et al., J. Phys. Chem. Ref. Data 4, 871-1178 (1975),
#      "Molten Salts: Volume 4, Part 2, Chlorides and Mixtures".
#      https://doi.org/10.1063/1.555527
# [R2] E. R. Van Artsdalen and I. S. Yaffe, J. Phys. Chem. 59, 118-127 (1955).
#      https://doi.org/10.1021/j150524a007
# [R3] P. Fellner, I. Votava, and M. Chrenkova-Paucirova, Chem. Zvesti 34,
#      330-334 (1980), "Density of the Molten System LiCl-NaCl-KCl".
#      External ternary validation/comparison source; not the implemented
#      default ternary equation.
#      https://www.chempap.org/index.php/chempap/article/view/1898
# [R4] F. Kohler, Monatsh. Chem. 91, 738-740 (1960).
#      https://doi.org/10.1007/BF00899814
# [R5] O. Redlich and A. T. Kister, Ind. Eng. Chem. 40, 345-348 (1948).
#      https://doi.org/10.1021/ie50458a036
# [R6] M. J. Mehl et al., Comput. Mater. Sci. 136, S1-S828 (2017).
#      https://doi.org/10.1016/j.commatsci.2017.01.017
#      Rock-salt prototype: https://aflow.org/p/AB_cF8_225_a_b-001/
# [R7] VASP Wiki, "POSCAR," https://vasp.at/wiki/POSCAR
# [R8] S. R. Hall, F. H. Allen, and I. D. Brown, Acta Cryst. A47, 655-685 (1991).
#      https://doi.org/10.1107/S010876739101067X
# [R9] LAMMPS documentation, https://docs.lammps.org/atom_style.html and
#      https://docs.lammps.org/read_data.html
# [R10] BIPM, SI Brochure, 9th ed.; N_A = 6.02214076e23 mol^-1 exactly.
#       https://doi.org/10.59161/AUEZ1291
# [R11] NIST Chemistry WebBook, SRD 69, formula/molecular-weight data.
#       https://doi.org/10.18434/T4D303
# [R12] G. J. Janz, "Thermodynamic and Transport Properties for Molten Salts:
#       Correlation Equations for Critically Evaluated Density, Surface Tension,
#       Electrical Conductance, and Viscosity Data," J. Phys. Chem. Ref. Data
#       17, Supplement 2 (1988). https://srd.nist.gov/JPCRD/jpcrdS2Vol17.pdf


if __name__ == "__main__":
    main()
