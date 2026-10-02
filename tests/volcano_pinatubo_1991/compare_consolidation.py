"""Compare a compiled pre-consolidation checkout with the current build.

Runs optical pre-runs and diagnostic fixtures in scratch space,
not the full formal driver or a background equilibrium simulation.
"""
import argparse
import filecmp
import itertools
import json
import math
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]


def run(command, work, timeout=180):
    result = subprocess.run(command, cwd=work, text=True, capture_output=True, timeout=timeout)
    if result.returncode:
        raise RuntimeError(result.stdout + result.stderr)
    return result


def inputs(root, work):
    for name in ('profile.dat', 'solar_flux.txt', 'xsecs', 'reactionsVerbatim.dat'):
        (work/name).symlink_to(root/'build'/name)


def compare(a, b):
    if filecmp.cmp(a, b, shallow=False):
        return {'byte_identical': True, 'max_relative_difference': 0.0}
    worst = 0.0
    with a.open() as first, b.open() as second:
        for lineno, (x, y) in enumerate(itertools.zip_longest(first, second), 1):
            if x == y:
                continue
            if x is None or y is None or x.startswith('#') or y.startswith('#'):
                raise AssertionError(f'{a.name}:{lineno}: changed header or row count')
            xx, yy = x.split(), y.split()
            if len(xx) != len(yy):
                raise AssertionError(f'{a.name}:{lineno}: changed column count')
            for lhs, rhs in zip(xx, yy):
                if lhs == rhs:
                    continue
                u, v = float(lhs), float(rhs)
                if not math.isfinite(u) or not math.isfinite(v):
                    raise AssertionError('Nonfinite comparison value')
                relative = abs(u-v)/max(abs(u), abs(v), 1e-99)
                worst = max(worst, relative)
                if relative > 5e-13:
                    raise AssertionError(f'{a.name}:{lineno}: relative difference {relative}')
    return {'byte_identical': False, 'max_relative_difference': worst}


def history_fixture(root, work, advance_chemistry):
    inputs(root, work)
    case = root/'tests/volcano_pinatubo_1991'
    original = (case/'test.f90').read_text()
    legacy = 'module volcano_history' in original
    program = '''program compare_history
  use patmo
  use patmo_volc
  use patmo_parameters
  use patmo_commons
  use patmo_constants
  USE_LEGACY
  implicit none
  real*8::t,dt,convergence
  integer::j,u
  call READ_SETTINGS()
  call patmo_init()
  call patmo_loadInitialProfile('profile.dat',unitH='km',unitX='1/cm3')
  call patmo_setFluxBB()
  call patmo_setGravity(9.8d2)
  wetdep=0d0
  va=0d0
  pa=0d0
  gd=0d0
  call patmo_volc_loadEvents('volcano_events.dat')
  ! Deliberately bypass spin-up only in this short equivalence fixture.
  call patmo_volc_startAfterEquilibrium()
  call HISTORY_CONFIGURE(200d0)
  call HISTORY_BEGIN()
  t=0d0
  convergence=100d0
  do while(t<117d0)
     dt=patmo_volc_limitStep(117d0-t)
     dt=HISTORY_LIMIT(t,dt)
     ADVANCE
     t=t+dt
     call HISTORY_SAMPLE(t)
  end do
  call HISTORY_CLOSE()
  open(newunit=u,file='fixture_final.dat',status='replace')
  do j=1,cellsNumber
     write(u,'(*(ES24.15E3,1X))') nall(j,:)
  end do
  close(u)
end program compare_history
'''
    old_names = ['read_run_settings', 'history_configure', 'history_begin',
                 'history_limit', 'history_sample', 'history_close']
    new_names = ['patmo_volc_readRunSettings', 'patmo_volc_historyConfigure',
                 'patmo_volc_historyBegin', 'patmo_volc_historyLimit',
                 'patmo_volc_historySample', 'patmo_volc_historyClose']
    keys = ['READ_SETTINGS', 'HISTORY_CONFIGURE', 'HISTORY_BEGIN',
            'HISTORY_LIMIT', 'HISTORY_SAMPLE', 'HISTORY_CLOSE']
    for key, name in zip(keys, old_names if legacy else new_names):
        program = program.replace(key, name)
    program = program.replace('USE_LEGACY',
                              'use volcano_history\n  use volcano_run_settings' if legacy else '')
    program = program.replace('ADVANCE', 'call patmo_run(dt,convergence)' if advance_chemistry
                              else 'call patmo_volc_setTime(t+dt)')
    (work/'check.f90').write_text(program)
    (work/'volcano_run.in').write_bytes((case/'volcano_run.in').read_bytes())
    (work/'volcano_events.dat').write_text(
        'start_s=17 duration_s=40 plume_center_km=22.5 plume_sigma_km=1.25 '
        'so2_column_cm2=1d15 ash_tau_550=0.4 ash_rise_time_s=30 '
        'sulfate_tau_550=0.15 sulfate_formation_day=25\n')
    (work/'volcano_history.in').write_text('''&history_output
 species_names='SO2','SO3', species_every_s=37, species_end_s=100,
 reaction_ids=0, reaction_every_s=41, reaction_end_s=100,
 solar_every_s=43, solar_end_s=100
/
''')
    names = ('opkda1 opkda2 opkdmain patmo_commons patmo_constants patmo_parameters '
             'patmo_utils patmo_rates patmo_reverseRates patmo_photo patmo_photoRates '
             'patmo_volc patmo_sparsity patmo_jacobian patmo_ode patmo').split()
    sources = []
    if legacy:
        (work/'history_modules.f90').write_text(original.split('\nprogram test\n', 1)[0])
        sources.append('history_modules.f90')
        names.insert(names.index('patmo_volc'), 'patmo_volc_transport')
    command = ['gfortran', '-O3', '-flto', '-funroll-loops', '-ffree-line-length-none',
               '-I'+str(root/'build')] + sources + ['check.f90']
    command += [str(root/'build'/(name+'.o')) for name in names] + ['-o', 'check']
    run(command, work)
    run([str(work/'check')], work)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-root', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--include-chemistry', action='store_true',
                        help='Also attempt a 117-second chemical fixture (180-second wall-time limit per run)')
    args = parser.parse_args()
    baseline = args.baseline_root.resolve()
    reference = (ROOT/'tests/volcano_pinatubo_1991/volcano_events.dat').read_text()
    cases = {
        'rise_hours': (reference, 'output_end_hour=12\noutput_time_step_hour=0.25\n'),
        'tail_days': (reference, 'output_end_day=5.25\noutput_time_step_day=0.5\n'),
        'auto_end': (reference, 'output_end_day=-1\noutput_time_step_day=2\ntau_floor=1d-6\n'),
        'delayed_seconds': (
            'start_s=17 duration_s=40 plume_center_km=22.5 plume_sigma_km=1.25 '
            'ash_tau_550=0.4 ash_rise_time_s=30 sulfate_tau_550=0.15 '
            'sulfate_start_delay_hour=0.01 sulfate_formation_day=1\n',
            'output_end_s=127\noutput_time_step_s=37\noutput_early_time_step_s=0\n'),
        'multiple_events': (reference + '\nstart_hour=12 duration_hour=0 '
            'plume_center_km=30 plume_sigma_km=1 ash_tau_550=0.1 ash_lifetime_day=1\n',
            'output_end_hour=48\noutput_time_step_hour=3\n'),
    }
    results = {}
    with tempfile.TemporaryDirectory(prefix='patmo_consolidation_compare_') as folder:
        scratch = Path(folder)
        for label, (event, config) in cases.items():
            paths = []
            for side, root in [('before', baseline), ('after', ROOT)]:
                work = scratch/(label+'_'+side)
                work.mkdir()
                inputs(root, work)
                (work/'volcano_events.dat').write_text(event)
                (work/'volcano_prerun.in').write_text('output_wavelength_step_nm=70\n'+config)
                run([str(root/'build/test_volcano')], work)
                paths.append(work)
            results[label] = {name: compare(paths[0]/name, paths[1]/name) for name in
                              ('volcano_optical_depth.dat', 'volcano_ash_profile.dat')}
            print(label + ': passed', flush=True)
        for chemistry in ((False, True) if args.include_chemistry else (False,)):
            label = 'short_chemistry' if chemistry else 'frozen_history'
            paths = []
            for side, root in [('before', baseline), ('after', ROOT)]:
                work = scratch/(label+'_'+side)
                work.mkdir()
                history_fixture(root, work, chemistry)
                paths.append(work)
            results[label] = {path.name: compare(path, paths[1]/path.name)
                              for path in paths[0].glob('volcano_history_*.dat')}
            if len(results[label]) != 5:
                raise AssertionError('Missing diagnostic output files')
            results[label]['fixture_final.dat'] = compare(paths[0]/'fixture_final.dat', paths[1]/'fixture_final.dat')
            print(label + ': passed', flush=True)
    args.report.write_text(json.dumps(results, indent=2)+'\n')
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
