"""Run after generation + make; never performs a full chemistry simulation."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CASE = Path(__file__).resolve().parent
BUILD = ROOT / 'build'
sys.path.insert(0, str(ROOT / 'src_py'))
from patmo_solar import HC_EV_NM, bin_average
import patmo_photoRates


class GenerationChecks(unittest.TestCase):
    def test_runtime_inputs_match_case(self):
        for name in (CASE / 'copylist.pcp').read_text().splitlines():
            name = name.strip()
            if name and not name.startswith('#'):
                self.assertEqual((CASE/name).read_bytes(), (BUILD/name).read_bytes(), name)

    def test_bin_average_linear_spectrum(self):
        low, high = np.array([250., 140.]), np.array([400., 250.])
        x = np.array([100., 150., 225., 400.])
        y = 2*x+3
        actual = bin_average(x, y, HC_EV_NM/high, HC_EV_NM/low)
        np.testing.assert_allclose(actual, low+high+3, rtol=1e-14)
        for badx, bady in (([100, 100, 400], [1, 2, 3]),
                           ([200, 400], [1, 2]), ([100, 400], [1, -1])):
            with self.assertRaises(ValueError):
                bin_average(badx, bady, HC_EV_NM/high, HC_EV_NM/low)

    def test_solar_order_widths_and_integrated_photons(self):
        metric = np.loadtxt(BUILD/'xsecs/photoMetric.dat')
        low, high = HC_EV_NM/metric[:,3], HC_EV_NM/metric[:,2]
        flux = np.loadtxt(BUILD/'solar_flux.txt')
        self.assertTrue(np.all(np.diff(high)<0))
        self.assertGreater(np.ptp(high-low), 0)
        source = pd.read_excel(CASE/'solar_flux.xlsx').iloc[:,:2].dropna().to_numpy(float)
        source = source[source[:,0].argsort()]
        x, y = source.T
        inside = x[(x>low.min()) & (x<high.max())]
        points = np.concatenate(([low.min()], inside, [high.max()]))
        values = np.interp(points, x, y)
        independent_integral = np.sum(np.diff(points)*0.5*(values[:-1]+values[1:]))
        self.assertAlmostEqual(np.dot(flux, high-low)/independent_integral, 1, places=12)
        expected = bin_average(x, y, metric[:,2], metric[:,3])
        np.testing.assert_allclose(flux, expected, rtol=2e-14)
        self.assertGreater(flux[0], flux[-1])  # This case's solar spectrum.

    def test_radiation_kernel(self):
        result = subprocess.run(['make', 'check_volcano_radiation'], cwd=BUILD,
                                text=True, capture_output=True, timeout=120)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)

    def test_radiation_settings_follow_generator(self):
        for angle in (45., 95.):
            with self.subTest(angle=angle), tempfile.TemporaryDirectory(prefix='patmo_angle_') as directory:
                work = Path(directory)
                generated = []
                def capture(source, destination, keys, values):
                    text = (ROOT/source).read_text()
                    for key, value in zip(keys, values):
                        text = text.replace(key, value)
                    generated.append(text)
                with patch.object(patmo_photoRates.patmo_string, 'fileReplaceBuild', side_effect=capture):
                    patmo_photoRates.buildPhotoRates(SimpleNamespace(photoReactions=[]),
                                                    SimpleNamespace(zenith_angle=angle, TOA_para=0.3))
                (work/'variant.f90').write_text(generated[0])
                expected = 3*np.exp(-0.2/np.cos(np.deg2rad(angle))) if angle < 90 else 0.
                literal = format(expected, '.16e').replace('e','d')
                (work/'check.f90').write_text(f'''program check
  use patmo_photoRates
  use patmo_parameters
  use patmo_commons
  implicit none
  real*8::tau(photoBinsNumber),flux(photoBinsNumber)
  energyLeft=10d0
  energyRight=20d0
  photoFlux=10d0
  tau=0.2d0
  flux=photoDirectFlux(tau)
  if(maxval(abs(flux-{literal}))>1d-12) error stop 'Wrong angle or TOA scaling'
end program check
''')
                command = ['gfortran','-O0','-fcheck=all','-I'+str(BUILD),
                           'variant.f90','check.f90']
                command += [str(BUILD/(n+'.o')) for n in ('patmo_commons','patmo_constants','patmo_parameters')]
                command += ['-o','check']
                result = subprocess.run(command,cwd=work,text=True,capture_output=True,timeout=60)
                self.assertEqual(result.returncode,0,result.stdout+result.stderr)
                result = subprocess.run([str(work/'check')],cwd=work,text=True,capture_output=True,timeout=60)
                self.assertEqual(result.returncode,0,result.stdout+result.stderr)

    def test_history_scheduler_and_spectral_closure(self):
        # Exercise the real diagnostic modules against a frozen test atmosphere.
        # Bypassing spin-up here is a unit-test fixture, not a science run.
        with tempfile.TemporaryDirectory(prefix='patmo_history_check_') as directory:
            work = Path(directory)
            (work/'check.f90').write_text('''program check_history
  use patmo
  use patmo_volc
  use patmo_parameters, only: krate
  implicit none
  real*8::t,dt
  call patmo_volc_readRunSettings()
  call patmo_init()
  call patmo_loadInitialProfile("profile.dat",unitH="km",unitX="1/cm3")
  call patmo_setFluxBB()
  call patmo_volc_loadEvents("volcano_events.dat")
  call patmo_volc_startAfterEquilibrium()
  call patmo_volc_historyConfigure(200d0)
  call patmo_volc_historyBegin()
  t=0d0
  do while(t<117d0)
     dt=patmo_volc_historyLimit(t,117d0-t)
     t=t+dt
     call patmo_volc_setTime(t)
     krate=0.123d0
     call patmo_volc_historySample(t)
     if(any(krate/=0.123d0)) error stop 'Diagnostics changed integrator rates'
  end do
  call patmo_volc_historyClose()
end program check_history
''')
            (work/'volcano_run.in').write_text((CASE/'volcano_run.in').read_text())
            (work/'volcano_events.dat').write_text(
                'start_s=17 duration_s=40 plume_center_km=22.5 plume_sigma_km=1.25 '
                'ash_tau_550=0.4 ash_rise_time_s=30\n')
            (work/'volcano_history.in').write_text('''&history_output
 species_names='SO2','SO3', species_every_s=37, species_end_s=100,
 reaction_ids=33,56, reaction_every_s=41, reaction_end_s=100,
 solar_every_s=43, solar_end_s=100
/
''')
            for name in ('profile.dat','solar_flux.txt','xsecs','reactionsVerbatim.dat'):
                (work/name).symlink_to(BUILD/name)
            names = ('opkda1 opkda2 opkdmain patmo_commons patmo_constants patmo_parameters '
                     'patmo_utils patmo_rates patmo_reverseRates patmo_photo patmo_photoRates '
                     'patmo_volc patmo_sparsity patmo_jacobian patmo_ode patmo').split()
            command = ['gfortran', '-O0', '-g', '-fcheck=all', '-ffree-line-length-none',
                       '-ffpe-trap=invalid,zero,overflow', '-I'+str(BUILD),
                       'check.f90']
            command += [str(BUILD/(name+'.o')) for name in names]
            command += ['-o', 'check']
            result = subprocess.run(command, cwd=work, text=True, capture_output=True, timeout=120)
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
            result = subprocess.run([str(work/'check')], cwd=work, text=True,
                                    capture_output=True, timeout=60)
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
            self.assertTrue((work/'volcano_history_species.dat').exists(), result.stdout+result.stderr)
            for suffix, expected in [('species',[0,37,74,100]), ('reactions',[0,41,82,100]),
                                     ('solar',[0,43,86,100])]:
                rows = np.loadtxt(work/f'volcano_history_{suffix}.dat', usecols=(0,))
                np.testing.assert_allclose(np.unique(rows), expected)
            wave = np.loadtxt(work/'volcano_history_wavelengths.dat')
            self.assertGreater(np.ptp(wave[:,3]), 0)
            self.assertTrue(np.all(wave[:,3]>0))
            metric = np.loadtxt(BUILD/'xsecs/photoMetric.dat')
            np.testing.assert_allclose(wave[:,1], HC_EV_NM/metric[:,0], rtol=1e-13)
            np.testing.assert_allclose(wave[:,3], HC_EV_NM*(1/metric[:,2]-1/metric[:,3]), rtol=1e-13)
            # Independent reconstruction from actual exported solar columns.
            solar = np.loadtxt(work/'volcano_history_solar.dat')
            reactions = np.loadtxt(work/'volcano_history_reactions.dat', dtype=str, skiprows=6)
            xsec = np.loadtxt(BUILD/'xsecs/SO2__SO_O.dat')[:,1]
            for time in (0,100):
                row = solar[(solar[:,0]==time) & (solar[:,1]==1)][0]
                reaction = reactions[(reactions[:,0].astype(float)==time) &
                                     (reactions[:,1]=='1') & (reactions[:,3]=='56')][0]
                self.assertEqual(reaction[4], 'J')
                expected = np.sum(row[3:]*wave[:,3]*xsec)
                self.assertAlmostEqual(float(reaction[5])/expected, 1, places=10)

    def test_wet_deposition_species_mapping(self):
        source = (CASE/'test.f90').read_text()
        setup = source.split('!calculate wet deposition\n',1)[1].split('  va(:) = 0d0',1)[0]
        routine = 'subroutine computewetdep' + source.split('\nsubroutine computewetdep',1)[1]
        routine = routine.split('end subroutine computewetdep',1)[0] + 'end subroutine computewetdep\n'
        with tempfile.TemporaryDirectory(prefix='patmo_wetdep_check_') as directory:
            work = Path(directory)
            for name in ('profile.dat','xsecs','reactionsVerbatim.dat'):
                (work/name).symlink_to(BUILD/name)
            # Execute the actual setup block from the case, not a duplicate list of calls.
            program = '''program check_wetdep
  use patmo
  use patmo_commons
  use patmo_parameters
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  logical::target(chemSpeciesNumber)
  real*8::saved(cellsNumber)
  integer::i
  call patmo_init()
  call patmo_loadInitialProfile("profile.dat",unitH="km",unitX="1/cm3")
  wetdep=0d0
''' + setup + '''
  target=.false.
  target([patmo_idx_COS,patmo_idx_CS2,patmo_idx_H2S,patmo_idx_SO2,patmo_idx_SO4,patmo_idx_H2SO4])=.true.
  if(.not.all(ieee_is_finite(wetdep))) error stop 'Nonfinite wet deposition'
  do i=1,chemSpeciesNumber
     if(target(i)) then
        if(any(wetdep(1:12,i)<=0d0)) error stop 'Missing intended wet deposition'
     else
        if(any(wetdep(:,i)/=0d0)) error stop 'Wet deposition on wrong species'
     end if
  end do
  if(any(wetdep(13:,:)/=0d0)) error stop 'Unexpected upper-layer wet deposition'
  if(abs(wetdep(1,patmo_idx_H2SO4)-dble(2.12E-05))>1d-15) error stop 'Changed H2SO4 value'
  ! COS kept its original index; also verify index-independent Henry-law numerics.
  saved=wetdep(:,patmo_idx_SO2)
  call computewetdep(patmo_idx_CS2,4.0d3)
  if(any(wetdep(:,patmo_idx_CS2)/=saved)) error stop 'Coefficient depends on species index'
  print *, 'Wet-deposition species mapping passed.'
end program check_wetdep
''' + routine
            (work/'check.f90').write_text(program)
            names = ('opkda1 opkda2 opkdmain patmo_commons patmo_constants patmo_parameters '
                     'patmo_utils patmo_rates patmo_reverseRates patmo_photo patmo_photoRates '
                     'patmo_volc patmo_sparsity patmo_jacobian patmo_ode patmo').split()
            command = ['gfortran','-O0','-g','-fcheck=all','-ffree-line-length-none',
                       '-ffpe-trap=invalid,zero,overflow','-I'+str(BUILD),'check.f90']
            command += [str(BUILD/(name+'.o')) for name in names]
            command += ['-o','check']
            result = subprocess.run(command,cwd=work,text=True,capture_output=True,timeout=120)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            result = subprocess.run([str(work/'check')],cwd=work,text=True,capture_output=True,timeout=60)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            self.assertIn('Wet-deposition species mapping passed.',result.stdout)

    def test_equilibrium_gate(self):
        with tempfile.TemporaryDirectory(prefix='patmo_equilibrium_check_') as directory:
            work = Path(directory)
            (work/'volcano_run.in').write_text('''&volcano_run
 spinup_max_years=1, volcano_duration_day=2, stable_days_required=3
/
''')
            (work/'volcano_events.dat').write_text(
                'start_s=0 duration_s=40 plume_center_km=22.5 plume_sigma_km=1.25 ash_tau_550=0.4\n')
            names = ('opkda1 opkda2 opkdmain patmo_commons patmo_constants patmo_parameters '
                     'patmo_utils patmo_rates patmo_reverseRates patmo_photo patmo_photoRates '
                     'patmo_volc patmo_sparsity patmo_jacobian patmo_ode patmo').split()
            command = ['gfortran', '-O0', '-g', '-fcheck=all', '-ffree-line-length-none',
                       '-ffpe-trap=invalid,zero,overflow', '-I'+str(BUILD),
                       str(CASE/'check_equilibrium.f90')]
            command += [str(BUILD/(name+'.o')) for name in names]
            command += ['-o', 'check']
            result = subprocess.run(command,cwd=work,text=True,capture_output=True,timeout=120)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            result = subprocess.run([str(work/'check')],cwd=work,text=True,capture_output=True,timeout=30)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            self.assertIn('Equilibrium gate checks passed.',result.stdout)
            result = subprocess.run([str(work/'check'),'fail'],cwd=work,text=True,capture_output=True,timeout=30)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('Volcano forcing was not started.',result.stdout)
            self.assertNotIn('Failure gate returned',result.stdout+result.stderr)


if __name__ == '__main__':
    unittest.main(verbosity=2)
