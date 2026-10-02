"""Regression checks after make: python3 tests/volcano_pinatubo_1991/test_prerun.py."""
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("volcano_viewer", ROOT / "tools/volcano_optical_depth_viewer.py")
viewer = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = viewer
spec.loader.exec_module(viewer)
viewer.configure_matplotlib(False)


class PreRunChecks(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="patmo_volcano_check_")
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)
        for name in ("profile.dat", "xsecs"):
            (self.work / name).symlink_to(ROOT / "build" / name)

    def run_case(self, event, end=48, step=1, early=0, floor=1e-6):
        (self.work / "volcano_events.dat").write_text(event + "\n")
        (self.work / "volcano_prerun.in").write_text(
            f"output_end_hour={end}\noutput_time_step_hour={step}\n"
            f"output_early_time_step_s={early}\noutput_early_until_hour=1\n"
            f"tau_floor={floor}\noutput_wavelength_step_nm=100\n"
        )
        result = subprocess.run([str(ROOT / "build/test_volcano")], cwd=self.work,
                                text=True, capture_output=True, timeout=90)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return viewer.load_volcano_file(self.work / "volcano_ash_profile.dat")

    def sulfate_case(self, extra=""):
        return ("start_day=0 duration_hour=0 plume_center_km=22.5 plume_sigma_km=1.25 "
                "sulfate_tau_550=1 sulfate_rise_time_hour=0 sulfate_formation_day=1 "
                "sulfate_lifetime_day=-1 sulfate_horizontal_lifetime_day=-1 "
                "sulfate_vertical_diffusion_cm2_s=0 " + extra)

    def field(self, data, name):
        return data.tau[:, :, data.field_names.index(name)]

    def test_formation_and_nonsettling_precursor(self):
        data = self.run_case(self.sulfate_case("sulfate_settling_cm_s=1"), end=24)
        p = self.field(data, "sulfate_precursor_potential")
        a = self.field(data, "sulfate_local_tau550")
        self.assertAlmostEqual(p[-1].sum(), np.exp(-1), places=8)
        self.assertAlmostEqual((p[-1] + a[-1]).sum(), 1, places=8)
        self.assertAlmostEqual(np.dot(p[-1], data.altitudes_km) / p[-1].sum(), 22.5, places=7)
        self.assertLess(np.dot(a[-1], data.altitudes_km) / a[-1].sum(), 22.2)
        np.testing.assert_allclose(self.field(data, "total_local_tau550"), a)
        self.assertNotIn("so2_source_cm-3_s-1", data.field_names)

    def test_uniform_source_analytic_solution(self):
        event = self.sulfate_case().replace("duration_hour=0", "duration_hour=9")
        data = self.run_case(event, end=24)
        expected = (24 / 9) * (np.exp(-15 / 24) - np.exp(-1))
        p = self.field(data, "sulfate_precursor_potential")[-1].sum()
        self.assertAlmostEqual(p, expected, delta=2e-7)
        self.assertAlmostEqual(self.field(data, "total_local_tau550")[-1].sum() + p, 1, places=8)

    def test_auto_end_includes_future_potential(self):
        event = self.sulfate_case().replace("sulfate_horizontal_lifetime_day=-1", "sulfate_horizontal_lifetime_day=2")
        data = self.run_case(event, end=-1, floor=1e-4)
        self.assertGreater(data.times[-1], 72)
        spectra = viewer.load_volcano_file(self.work / "volcano_optical_depth.dat")
        self.assertLess(np.max(spectra.tau[-1]), 1e-4)
        self.assertGreater(np.max(spectra.tau[1:]), 1e-4)

    def test_early_rise_and_complete_spectra(self):
        event = Path(__file__).with_name("volcano_events.dat").read_text()
        data = self.run_case(event, end=96, early=300)
        self.assertEqual(data.times[-1], 96)
        self.assertTrue(np.any((data.times > 0) & (data.times < 0.5)))
        ash = self.field(data, "ash_local_tau550")
        means = ash @ data.altitudes_km / np.maximum(ash.sum(axis=1), 1e-99)
        self.assertLess(means[1], 10)
        # Continuous supply leaves material throughout the column, so its mean
        # height after one transit is not the leading-edge height.
        self.assertGreater(means[np.argmin(abs(data.times - 0.5))], 8)
        self.assertGreater(ash[np.argmin(abs(data.times - 0.5)), data.altitudes_km >= 20].sum(), 1e-5)
        self.assertGreater(ash[np.argmin(abs(data.times - 9))].sum(), ash[1].sum())
        self.assertGreater(np.sum(self.field(data, "sulfate_local_tau550")[-1]), 1e-6)
        self.assertTrue(np.all(data.tau >= 0))
        spectra = viewer.load_volcano_file(self.work / "volcano_optical_depth.dat")
        np.testing.assert_array_equal(spectra.times, data.times)
        # Output uses eight significant digits; account for independent rounding.
        np.testing.assert_allclose(spectra.tau[:, :, 0], spectra.tau[:, :, 1] + spectra.tau[:, :, 2], rtol=1e-7, atol=1e-12)
        self.assertTrue(np.all(np.isfinite(spectra.tau[-1, :, :, 0])))
        options = viewer.PlotOptions(field_name="sulfate_local_tau550", time_min=0, time_max=96)
        viewer.save_heatmap(data, options, self.work / "ssa.png")
        self.assertGreater((self.work / "ssa.png").stat().st_size, 10000)
        with self.assertRaises(ValueError):
            viewer.time_range_indices(data, viewer.PlotOptions(time_max=100))

    def test_output_sampling_does_not_change_state(self):
        event = self.sulfate_case("sulfate_settling_cm_s=1")
        hourly = self.run_case(event, end=24, step=1)
        half = self.run_case(event, end=24, step=0.5)
        np.testing.assert_allclose(hourly.tau[-1], half.tau[-1], rtol=3e-4, atol=1e-7)

    def rising_case(self):
        return ("start_hour=0 duration_hour=9 plume_bottom_km=20 plume_top_km=25 "
                "ash_tau_550=1 ash_rise_start_km=0 ash_rise_time_hour=0.5 "
                "ash_lifetime_day=-1 ash_horizontal_lifetime_day=-1 "
                "ash_vertical_diffusion_cm2_s=0 ash_settling_cm_s=0 "
                "sulfate_tau_550=1 sulfate_rise_start_km=0 sulfate_rise_time_hour=0.5 "
                "sulfate_formation_day=1 sulfate_lifetime_day=-1 "
                "sulfate_horizontal_lifetime_day=-1 sulfate_vertical_diffusion_cm2_s=0 "
                "sulfate_settling_cm_s=0")

    def test_sustained_vent_source_and_bounded_ascent(self):
        data = self.run_case(self.rising_case(), end=12, step=0.25)
        ash = self.field(data, "ash_local_tau550")
        p = self.field(data, "sulfate_precursor_potential")
        a = self.field(data, "sulfate_local_tau550")
        expected = np.minimum(data.times / 9, 1)
        np.testing.assert_allclose(ash.sum(axis=1), expected, atol=2e-8)
        np.testing.assert_allclose((p+a).sum(axis=1), expected, atol=2e-8)
        # The vent and intermediate layers are still fed near the end, not just
        # during the initial half hour. All upward transport stops at the band.
        late = np.argmin(abs(data.times - 8))
        self.assertTrue(np.all(ash[late, data.altitudes_km <= 20] > 1e-5))
        self.assertTrue(np.all(p[late, data.altitudes_km <= 20] > 1e-5))
        self.assertLess(np.max(ash[:, data.altitudes_km > 25]), 1e-14)
        self.assertLess(np.max(p[:, data.altitudes_km > 25]), 1e-14)
        self.assertLess(np.max(a[:, data.altitudes_km > 25]), 1e-14)
        self.assertLess(ash[-1, data.altitudes_km <= 10].sum(), 1e-5)
        # After the clearing interval, ash with zero settling/mixing is fixed.
        np.testing.assert_allclose(ash[-1], ash[np.argmin(abs(data.times-10))], atol=1e-12)
        self.assertLess(a[late].sum(), ash[late].sum())
        expected_p = (24/9) * (np.exp(-(12-9)/24) - np.exp(-12/24))
        self.assertAlmostEqual(p[-1].sum(), expected_p, delta=2e-7)

    def test_delayed_continuous_source_and_sampling(self):
        event = self.rising_case().replace("start_hour=0", "start_hour=2")
        event += " sulfate_start_delay_hour=1"
        data = self.run_case(event, end=12, step=0.25)
        p = self.field(data, "sulfate_precursor_potential")
        a = self.field(data, "sulfate_local_tau550")
        expected = np.clip((data.times-1)/9, 0, 1)
        np.testing.assert_allclose((p+a).sum(axis=1), expected, atol=2e-8)
        coarse = self.run_case(event, end=12, step=1)
        np.testing.assert_allclose(data.tau[-1], coarse.tau[-1], rtol=3e-4, atol=1e-7)

    def test_formed_ssa_settles_independently_after_ascent(self):
        # Immediate conversion is a limiting test, not the reference scenario.
        event = self.rising_case().replace("sulfate_formation_day=1", "sulfate_formation_day=0")
        event = event.replace("sulfate_settling_cm_s=0", "sulfate_settling_cm_s=1")
        data = self.run_case(event, end=24, step=1)
        ash = self.field(data, "ash_local_tau550")
        a = self.field(data, "sulfate_local_tau550")
        self.assertEqual(self.field(data, "sulfate_precursor_potential").max(), 0)
        mean = a @ data.altitudes_km / np.maximum(a.sum(axis=1), 1e-99)
        at10 = np.argmin(abs(data.times-10))
        # A sharp, one-cell accumulation has sub-cell reconstruction error;
        # the transport kernel separately tests a smooth profile analytically.
        self.assertAlmostEqual(mean[-1]-mean[at10], -1*14*3600/1e5, delta=0.1)
        np.testing.assert_allclose(ash[-1], ash[at10], atol=1e-12)
        self.assertLess(a[-1].sum(), 1)  # Physical vent-layer settling outflow.
        self.assertGreater(a[-1].sum(), 0.99)

    def test_top_layer_is_in_shielding(self):
        event = ("start_day=0 duration_hour=0 plume_center_km=60 plume_sigma_km=0.01 "
                 "ash_tau_550=1 ash_rise_time_hour=0 ash_lifetime_day=-1 "
                 "ash_horizontal_lifetime_day=-1 ash_vertical_diffusion_cm2_s=0")
        data = self.run_case(event, end=1)
        local = self.field(data, "ash_local_tau550")[-1]
        cumulative = self.field(data, "ash_column_tau550")[-1]
        self.assertAlmostEqual(local[-1], 1)
        np.testing.assert_allclose(cumulative, np.cumsum(local[::-1])[::-1])

    def test_display_threshold_before_normalization(self):
        data = self.run_case(self.sulfate_case(), end=1)
        values = np.zeros((2, len(data.layers)))
        values[0, 20] = 1e-14
        values[1, 20] = 0.1
        plotted, _ = viewer.apply_vertical_view(data, values, "SSA", viewer.PlotOptions(normalize_profile=True))
        self.assertTrue(np.all(np.isnan(plotted[0])))
        self.assertAlmostEqual(np.nansum(plotted[1]), 1)
        plotted, _ = viewer.apply_vertical_view(data, values, "SSA", viewer.PlotOptions(column_floor=0))
        self.assertEqual(plotted[0, 20], 1e-14)
        viewer.save_heatmap(data, viewer.PlotOptions(column_floor=10), self.work / "empty.png")
        viewer.make_animation(data, viewer.PlotOptions(dpi=60), self.work / "column.gif")
        self.assertGreater((self.work / "column.gif").stat().st_size, 1000)

    def test_reject_incomplete_and_duplicate_data(self):
        path = self.work / "bad.dat"
        header = "# columns: time layer altitude_km ash_local_tau550\n"
        for rows in ("0 1 1 0\n0 2 2 0\n1 1 1 0\n",
                     "0 1 1 0\n0 1 1 1\n", "0 1 1 0\n1 1 1\n"):
            path.write_text(header + rows)
            with self.assertRaises(ValueError):
                viewer.load_volcano_file(path)


if __name__ == "__main__":
    unittest.main(verbosity=2)
