import importlib.util
import math
import random
import tempfile
import unittest
import sys
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "build_mixed_alkali_chloride_rocksalt.py"
spec = importlib.util.spec_from_file_location("builder", SCRIPT)
builder = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = builder
assert spec.loader is not None
spec.loader.exec_module(builder)


class BuilderTests(unittest.TestCase):
    def test_legacy_li_k_volume_is_unchanged(self):
        vm = builder.mixture_molar_volume_cm3_mol(
            1073.15,
            63 / 108,
            0.0,
            45 / 108,
            use_excess=True,
        )
        self.assertAlmostEqual(vm, 37.6059422756, places=9)


    def test_legacy_linak_ternary_volume_field_is_frozen(self):
        # Compatibility oracle only: these values freeze released v1.3.0 behavior.
        # They are not empirical-validity oracles and must not be interpreted as
        # reproducing the distinct Fellner et al. ternary equation.
        T = 1073.15
        anchors = [
            ((1 / 6, 1 / 6, 2 / 3), 45.273151240704706),
            ((1 / 6, 2 / 3, 1 / 6), 38.313467807963427),
            ((1 / 3, 1 / 3, 1 / 3), 39.018773991093049),
            ((2 / 3, 1 / 6, 1 / 6), 33.677252249816618),
            ((63 / 108, 0.0, 45 / 108), 37.605942275646626),
        ]
        for (x_li, x_na, x_k), expected_vm in anchors:
            got = builder.mixture_molar_volume_cm3_mol(
                T, x_li, x_na, x_k, use_excess=True
            )
            self.assertAlmostEqual(got, expected_vm, places=12)

    def test_legacy_pure_density_independent_anchors(self):
        # Independent source-table anchors; expected values are not generated
        # from production coefficients.
        anchors = [
            ("LiCl", 1000.0, 1.4478000000000000),
            ("NaCl", 1150.0, 1.5129400000000000),
            ("KCl", 1080.0, 1.5069701201600000),
        ]
        for salt, T_K, expected in anchors:
            self.assertAlmostEqual(builder.rho_pure_salt(T_K, salt), expected, places=12, msg=salt)

    def test_extended_pure_density_independent_anchors(self):
        # Janz 1988 values transcribed independently from the source table.
        anchors = {
            "RbCl": (715.0, 2.248, 0.000883),
            "CsCl": (645.0, 2.790, 0.001065),
            "BeCl2": (415.0, 1.540, 0.001100),
            "MgCl2": (714.0, 1.680, 0.000271),
            "CaCl2": (775.0, 2.085, 0.000422),
            "SrCl2": (874.0, 2.727, 0.000578),
            "BaCl2": (962.0, 3.174, 0.000681),
        }
        for salt, (tm_c, rho_m, k) in anchors.items():
            self.assertAlmostEqual(builder.rho_pure_salt(tm_c + 273.15, salt), rho_m, places=12, msg=salt)
            self.assertAlmostEqual(
                builder.rho_pure_salt(tm_c + 373.15, salt),
                rho_m - 100.0 * k,
                places=12,
                msg=f"{salt} slope",
            )

    def test_legacy_pure_temperature_range_warnings(self):
        pure_cases = [
            ("Li", 880.0, 1070.0, 1200.0),
            ("Na", 1090.0, 1220.0, 900.0),
            ("K", 1060.0, 1200.0, 900.0),
        ]
        for species, lo, hi, outside in pure_cases:
            fractions = {s: 0.0 for s in builder.CATION_ORDER}
            fractions[species] = 1.0
            self.assertEqual(builder.warn_temperature_ranges(lo, fractions), [], msg=f"{species} low endpoint")
            self.assertEqual(builder.warn_temperature_ranges(hi, fractions), [], msg=f"{species} high endpoint")
            warnings = builder.warn_temperature_ranges(outside, fractions)
            self.assertTrue(any(f"pure-{builder.SALT_FORMULA[species]}" in w for w in warnings), msg=species)


    def test_mixture_checks_active_legacy_pure_component_ranges(self):
        fractions = {s: 0.0 for s in builder.CATION_ORDER}
        fractions.update({"Li": 0.25, "Cs": 0.25, "Mg": 0.25, "Ca": 0.25})
        warnings = builder.warn_temperature_ranges(1090.0, fractions)
        self.assertTrue(any("pure-LiCl" in w for w in warnings))

    def test_lina_anchor_temperature_range_warning(self):
        fractions = {s: 0.0 for s in builder.CATION_ORDER}
        fractions.update({"Li": 0.5, "Na": 0.5})
        at_endpoint = builder.warn_temperature_ranges(1123.0, fractions)
        self.assertFalse(any("50:50 source-correlation range" in w for w in at_endpoint))
        warnings = builder.warn_temperature_ranges(1150.0, fractions)
        self.assertTrue(any("50:50 source-correlation range" in w for w in warnings))

    def test_ideal_only_suppresses_inactive_legacy_pair_warning(self):
        fractions = {s: 0.0 for s in builder.CATION_ORDER}
        fractions.update({"Li": 0.5, "Na": 0.5})
        warnings = builder.warn_temperature_ranges(850.0, fractions, use_excess=False)
        self.assertTrue(any("pure-LiCl" in w for w in warnings))
        self.assertTrue(any("pure-NaCl" in w for w in warnings))
        self.assertFalse(any("legacy binary model" in w for w in warnings))

    def test_cli_ideal_only_reports_only_active_validity_models(self):
        with tempfile.TemporaryDirectory() as td:
            prefix = Path(td) / "ideal"
            proc = subprocess.run(
                [
                    sys.executable, str(SCRIPT),
                    "--nx", "1", "--ny", "1", "--nz", "1",
                    "--temp-k", "850",
                    "--x-li", "0.5", "--x-na", "0.5", "--x-k", "0",
                    "--ideal-only",
                    "--format", "lammps",
                    "--output-prefix", str(prefix),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
        self.assertIn("Excess-volume terms", proc.stdout)
        self.assertIn(": off, pure-component molar-volume interpolation only", proc.stdout)
        self.assertIn("pure-LiCl", proc.stdout)
        self.assertIn("pure-NaCl", proc.stdout)
        self.assertNotIn("legacy binary model", proc.stdout)

    def test_cli_active_legacy_binary_warnings_are_not_over_suppressed(self):
        # Real-owner positive counterexamples for the assembled CLI path.
        # These specifically guard main() propagation: helper-only tests would
        # remain green if the CLI accidentally forced use_excess=False for
        # diagnostics while still evaluating excess terms for the volume.
        pair_cases = [
            ("LiCl-NaCl", 850.0, {"li": 0.5, "na": 0.5, "k": 0.0}),
            ("NaCl-KCl", 900.0, {"li": 0.0, "na": 0.5, "k": 0.5}),
            ("KCl-LiCl", 650.0, {"li": 0.5, "na": 0.0, "k": 0.5}),
        ]
        for pair_name, temp_k, fractions in pair_cases:
            with self.subTest(pair=pair_name), tempfile.TemporaryDirectory() as td:
                prefix = Path(td) / "active"
                proc = subprocess.run(
                    [
                        sys.executable, str(SCRIPT),
                        "--nx", "1", "--ny", "1", "--nz", "1",
                        "--temp-k", str(temp_k),
                        "--x-li", str(fractions["li"]),
                        "--x-na", str(fractions["na"]),
                        "--x-k", str(fractions["k"]),
                        "--format", "lammps",
                        "--output-prefix", str(prefix),
                    ],
                    check=True,
                    capture_output=True,
                    text=True,
                )
            self.assertIn("Excess-volume terms", proc.stdout)
            self.assertNotIn(": off, pure-component molar-volume interpolation only", proc.stdout)
            self.assertIn(pair_name, proc.stdout)
            self.assertIn("legacy binary model", proc.stdout)

    def test_temperature_range_endpoint_roundoff_is_tolerated(self):
        fractions = {s: 0.0 for s in builder.CATION_ORDER}
        fractions["Mg"] = 1.0
        # 826 C expressed through Kelvin used to spuriously compare just above
        # the tabulated upper endpoint after binary floating-point subtraction.
        self.assertEqual(builder.warn_temperature_ranges(826.0 + 273.15, fractions), [])
        warnings = builder.warn_temperature_ranges(826.0 + 273.15 + 1.0e-6, fractions)
        self.assertTrue(any("pure-MgCl2" in w for w in warnings))

    def test_ideal_volume_is_composition_interpolation(self):
        T = 1090.0
        T_C = T - 273.15
        fractions = {s: 0.0 for s in builder.CATION_ORDER}
        fractions["Rb"] = 0.5
        fractions["Mg"] = 0.5

        # Independent oracle from source coefficients and formula masses; do
        # not call production pure-volume functions to build the expected value.
        rho_rbcl = 2.248 - 0.000883 * (T_C - 715.0)
        rho_mgcl2 = 1.680 - 0.000271 * (T_C - 714.0)
        m_rbcl = 85.46780 + 35.45300
        m_mgcl2 = 24.30500 + 2.0 * 35.45300
        expected = 0.5 * (m_rbcl / rho_rbcl) + 0.5 * (m_mgcl2 / rho_mgcl2)

        got = builder.mixture_molar_volume_from_fractions(T, fractions, use_excess=False)
        self.assertAlmostEqual(got, expected, places=12)

    def test_charge_balance_mixed_monovalent_divalent(self):
        counts = {s: 0 for s in builder.CATION_ORDER}
        counts.update({"Li": 1, "Cs": 1, "Mg": 1, "Ca": 1})
        self.assertEqual(builder.required_chloride_count(counts), 6)
        cations = ["Li", "Cs", "Mg", "Ca"]
        atoms = builder.build_mixed_chloride_atoms(1, 1, 1, cations, random.Random(12345))
        self.assertEqual(sum(a.symbol == "Cl" for a in atoms), 6)
        self.assertEqual(len(atoms), 10)
        q = sum(builder.CHARGES[a.symbol] for a in atoms)
        self.assertAlmostEqual(q, 0.0, places=12)

    def test_pure_divalent_has_two_chlorides_per_cation(self):
        cations = ["Ca"] * 4
        atoms = builder.build_mixed_chloride_atoms(1, 1, 1, cations, random.Random(7))
        self.assertEqual(sum(a.symbol == "Ca" for a in atoms), 4)
        self.assertEqual(sum(a.symbol == "Cl" for a in atoms), 8)
        self.assertEqual(len(atoms), 12)

    def test_fraction_rounding_preserves_cation_site_count(self):
        fractions = {s: 0.0 for s in builder.CATION_ORDER}
        fractions.update({"Li": 0.17, "Rb": 0.19, "Mg": 0.31, "Ba": 0.33})
        counts = builder.counts_from_fraction_map(108, fractions)
        self.assertEqual(sum(counts.values()), 108)
        self.assertTrue(all(v >= 0 for v in counts.values()))

    def test_lammps_charge_output_contains_plus_two_and_balances(self):
        counts = {s: 0 for s in builder.CATION_ORDER}
        counts.update({"Na": 2, "Sr": 2})
        lattice = builder.lattice_vectors(10.0, 1, 1, 1)
        atoms = builder.build_mixed_chloride_atoms(1, 1, 1, ["Na", "Na", "Sr", "Sr"], random.Random(2))
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / "cell.data"
            builder.write_lammps_data(out, atoms, lattice, "test", "charge")
            text = out.read_text()
        self.assertIn("2.000000", text)
        self.assertEqual(sum(a.symbol == "Cl" for a in atoms), 6)
        self.assertAlmostEqual(sum(builder.CHARGES[a.symbol] for a in atoms), 0.0, places=12)

    def test_all_extended_pure_range_endpoints_and_outside(self):
        # Independent source ranges, expressed in Celsius.
        ranges_c = {
            "Rb": (715.0, 923.0),
            "Cs": (645.0, 906.0),
            "Be": (415.0, 473.0),
            "Mg": (714.0, 826.0),
            "Ca": (775.0, 950.0),
            "Sr": (874.0, 1037.0),
            "Ba": (962.0, 1081.0),
        }
        for species, (lo_c, hi_c) in ranges_c.items():
            fractions = {s: 0.0 for s in builder.CATION_ORDER}
            fractions[species] = 1.0
            self.assertEqual(builder.warn_temperature_ranges(lo_c + 273.15, fractions), [], msg=f"{species} low")
            self.assertEqual(builder.warn_temperature_ranges(hi_c + 273.15, fractions), [], msg=f"{species} high")
            below = builder.warn_temperature_ranges(lo_c + 273.15 - 1.0e-6, fractions)
            above = builder.warn_temperature_ranges(hi_c + 273.15 + 1.0e-6, fractions)
            salt = builder.SALT_FORMULA[species]
            self.assertTrue(any(f"pure-{salt}" in w for w in below), msg=f"{species} below")
            self.assertTrue(any(f"pure-{salt}" in w for w in above), msg=f"{species} above")

    def test_legacy_pair_range_endpoints_and_outside(self):
        pair_cases = [
            (("Li", "Na"), "LiCl-NaCl", 893.0, 1123.0),
            (("Na", "K"), "NaCl-KCl", 943.0, 1300.0),
            (("K", "Li"), "KCl-LiCl", 668.0, 1212.0),
        ]
        for (a, b), pair_name, lo, hi in pair_cases:
            fractions = {s: 0.0 for s in builder.CATION_ORDER}
            fractions[a] = 0.5
            fractions[b] = 0.5
            at_lo = builder.warn_temperature_ranges(lo, fractions)
            at_hi = builder.warn_temperature_ranges(hi, fractions)
            self.assertFalse(any(pair_name in w and "legacy binary model" in w for w in at_lo), msg=f"{pair_name} low")
            self.assertFalse(any(pair_name in w and "legacy binary model" in w for w in at_hi), msg=f"{pair_name} high")
            below = builder.warn_temperature_ranges(lo - 1.0e-6, fractions)
            above = builder.warn_temperature_ranges(hi + 1.0e-6, fractions)
            self.assertTrue(any(pair_name in w and "legacy binary model" in w for w in below), msg=f"{pair_name} below")
            self.assertTrue(any(pair_name in w and "legacy binary model" in w for w in above), msg=f"{pair_name} above")

    def test_cli_finite_cell_volume_independent_anchor(self):
        with tempfile.TemporaryDirectory() as td:
            prefix = Path(td) / "legacy"
            proc = subprocess.run(
                [
                    sys.executable, str(SCRIPT),
                    "--nx", "3", "--ny", "3", "--nz", "3",
                    "--temp-c", "800",
                    "--x-li", "0.585", "--x-k", "0.415",
                    "--format", "lammps",
                    "--output-prefix", str(prefix),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
        line = next(line for line in proc.stdout.splitlines() if line.startswith("Equilibrium cell volume"))
        value = float(line.split(":", 1)[1].split()[0])
        self.assertAlmostEqual(value, 6744.1827211123, places=9)

    def test_legacy_monovalent_site_geometry_is_identical(self):
        cations = ["Li", "Na", "K", "Li"]
        old = builder.build_rocksalt_atoms(1, 1, 1, cations)
        new = builder.build_mixed_chloride_atoms(1, 1, 1, cations, random.Random(99))
        self.assertEqual([(a.symbol, a.frac) for a in old], [(a.symbol, a.frac) for a in new])


if __name__ == "__main__":
    unittest.main()
