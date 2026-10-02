# Pinatubo 1991 Column Test

This Pinatubo-informed sensitivity experiment uses the chemistry network and
background profile originally derived from modern_sulfur_cycle. It is not an observation-fitted
Pinatubo reconstruction. SSA here means sulfate aerosol, not single-scattering
albedo.

## First-stage Physics

Three layer-integrated budgets are tracked independently: ash optical depth,
unconverted sulfate precursor potential, and formed SSA optical depth.
Precursor potential is a prescribed tau550-equivalent quantity, not gas SO2
concentration. It does not enter radiative shielding and does not sediment.
It rises and mixes with the air, and is converted by a first-order law:

```text
precursor(t+dt) = precursor(t) * exp(-dt / formation_time)
formed_SSA     = precursor(t) - precursor(t+dt)
total_tau      = ash_tau + formed_SSA_tau
```

Only formed SSA undergoes SSA settling and additional particle removal.
Horizontal export acts on both precursor and formed SSA with the same
sulfate horizontal lifetime. All three budgets use finite-volume transport:
MC-limited advection, SSP-RK2 integration and density-weighted mixing flux

```text
flux_mix = -Kzz * n_air * d[(layer_tau / dz) / n_air] / dz
```

Each internal face has one flux, shared by its two neighboring layers.
Downward bottom flux leaves the column. The top has zero flux. Internal
substeps enforce stability. No arbitrary removal is imposed at high altitude.
These are optical-budget conservation statements, not sulfur-mass closure.
Fixed extinction efficiency is an assumption.

For positive rise_time, ash and sulfate precursor enter a fixed vent layer
at a constant rate throughout the single event. A vent below the grid maps
to the lowest layer. No source moves to the stratosphere after startup.
The prescribed upward speed is (center - rise_start) / rise_time below the
target band and decreases linearly to zero between center - 2 sigma and
center + 2 sigma. The taper's lower height cannot be below the vent.
Ascent operates throughout emission plus one nominal transit time to carry
the final emissions upward; it then switches off. This extra interval adds
no source. The taper means rise_time is not an exact arrival time.

Precursor has no gravitational settling. SSA forms in its current layer;
formed particles share the sulfate plume airflow, with their own settling
and mixing. No ash profile is copied into SSA, and no direct SSA source is
added. Ash and sulfate retain their separate transport parameters. After
the ascent interval, only mixing, settling, formation and losses remain.
Particles can mix above the target band, but ascent cannot pump them there.
Zero rise_time explicitly retains Gaussian high-altitude injection.

This is an effective 1-D sustained-column prescription, not a buoyant plume
solver: it does not resolve entrainment, thermal energy, overshoot or lateral
detrainment. Re-entrainment within the active column is implicit. The target
band now defines the ascent taper, not a guaranteed Gaussian particle cloud.
The chemical SO2 source is unchanged (Gaussian high-altitude injection);
the optical precursor remains independent of chemical SO2 and sulfur mass.

## Current Reference Settings

| Input | Default | Interpretation |
| --- | --- | --- |
| duration_hour | 9 | Averaged climactic injection duration |
| ash/sulfate_rise_time_hour | 0.5 | Nominal ascent timescale, independent of duration |
| plume_bottom/top_km | 20 / 25 | Ascent taper band; Gaussian SO2 injection interval |
| ash/sulfate_rise_start_km | 0 | Idealized ground coordinate, not measured vent elevation |
| so2_mass_tg / plume_radius_km | 20 / 2500 | Chemistry source over a large regional averaging area |
| ash_tau_550 / sulfate_tau_550 | 0.4 / 0.15 | Prescribed optical budgets, not realized local AOD peaks |
| ash/sulfate particle radius_um | 2 / 0.4 | Fixed spherical fine-ash / mature-SSA sensitivity radii |
| ash/sulfate Cc | 1.9 / 6.3 | Fixed reference-height slip factors, not altitude-dependent |
| ash/sulfate settling_km_day | 0.21044 / 0.017863 | Derived from respective radius, density, viscosity and Cc |
| ash/sulfate vertical_diffusion_cm2_s | background | Same air Kzz profile; numeric overrides remain available |
| sulfate_formation_day | 25 | Conversion e-folding time, not completion time |
| ash/sulfate horizontal_lifetime_day | 1 / 14 | Local export assumptions, not global residence times |
| ash_lifetime_day | 2 | Additional residual ash loss after emission ends |
| sulfate_lifetime_day | -1 | No additional residual loss on top of explicit processes |

The 9-hour duration and 25-day timescale are informed by
[Guo et al. (2004), introduction, abstract and Table 5](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2003GC000654).
The observed SO2 removal e-folding time is used as an approximate conversion
parameter; it is not a uniquely measured chemical rate. The assumed 0.5-hour
rise should be varied over 0.25-1 hour as a sensitivity test.

The radius 2500 km defines regional averaging, not a near-vent location.
Local export times and optical budgets require observations from the same
spatial footprint before quantitative comparison. Fixed Cc, viscosity and
radius are reference-height approximations. Altitude-dependent drag, particle
size evolution, ice/ash aggregation, sulfur mass coupling and multiple
scattering are not implemented. A small nonzero top-layer tail is not proof
of a significant high-altitude aerosol layer.

See reports/pinatubo_parameter_review_2026-09-08.md for the detailed prior
audit and primary references. Its earlier algorithm descriptions are historical;
this README describes the revised implementation.

## Source Layout

The authoritative volcanic implementation is src_f90/patmo_volc.f90. It owns
event parsing, run settings, the equilibrium monitor, ash/precursor/SSA
transport and formation, optical shielding, standalone pre-run configuration,
and formal species/reaction/solar diagnostics. The transport routines are
part of this module, not a separate patmo_volc_transport compilation unit.

The case test.f90 retains background initialization, wet deposition, the
PATMO integration loops and general final outputs. It calls the volcanic
module to check equilibrium, control forcing and sample diagnostics.
test_volcano.f90 only configures the pre-run, loads the profile/optical grid
and calls the volcanic module. The module does not depend on the high-level
patmo module, so there is no circular compilation dependency.

All four volcanic input files, their defaults, physical equations, output
filenames, columns and units retain their pre-consolidation behavior.
The desktop viewer does not need changes for this structural consolidation.

## Build and Pre-run

From the repository root, regenerate inputs/code, then compile:

```bash
./compile.sh volcano_pinatubo_1991
cd build
make
./test_volcano
```

For conversion, generation and compilation in one command, run:

```bash
./tests/volcano_pinatubo_1991/compile_volcano_pinatubo_1991.sh
```

This case-local script works from any working directory. It delegates to the
same root script, copies the case test.f90 verbatim, verifies all manifest
inputs against build, then runs make. It stops on failure and never starts
test or test_volcano. Neither entry rewrites the case copylist.pcp or test.f90.
The root compile.sh still only generates the model, without running make. A direct
python3 patmo -test=volcano_pinatubo_1991 uses existing options.opt,
reaction_network.ntw and profile.dat; compile.sh first refreshes them from
their Excel sources. Both generator routes align solar_flux.xlsx to the
actual generated photoMetric.dat, storing photons cm-2 s-1 nm-1 as bin averages.
Do not edit build as the only copy of a lasting change.

## Input Ownership

| File | Purpose | Consumed by |
| --- | --- | --- |
| settings.xlsx | Grid, radiation angle/scaling, chemical boundary settings | compile.sh -> options.opt |
| reaction_network.xlsx | Chemical network | compile.sh -> reaction_network.ntw |
| profile.xlsx | Initial background profile | compile.sh -> profile.dat |
| solar_flux.xlsx | Physical wavelength/photon spectral flux | generator -> solar_flux.txt |
| volcano_events.dat | Single averaged event, ash and SSA parameters | Both executables |
| volcano_prerun.in | Pre-run time/wavelength sampling and end threshold | test_volcano |
| volcano_run.in | Spin-up upper limit, equilibrium criteria, model duration | test |
| volcano_history.in | Selected species/reactions and output intervals/windows | test |
| copylist.pcp | Authoritative list of runtime inputs and drivers | generator |

partial_H2SO4.txt, vapor_H2SO4.txt and SO4_deposition_rate.txt remain
physical input tables. Rainout-*.txt diagnostics are not copied into build;
the driver generates the active rainout diagnostics. Existing copies in this
case directory are historical references, not authoritative inputs.

Solar radiation uses actual per-bin wavelength widths, not one constant
delta-lambda. The photochemistry and history output share photoDirectFlux
and photoWavelengthWidths. Zenith angles >= 90 degrees have zero direct
flux; diffuse/multiple scattering is not modeled. Wet deposition uses generated
patmo_idx_* names rather than legacy numeric species positions. Changing the
reaction network still requires choosing appropriate Henry constants and
isotopologue deposition assumptions; those are not inferred automatically.

The standalone pre-run does not call the chemical ODE. It reads
profile.dat, xsecs/photoMetric.dat and volcano_events.dat. Settings in
volcano_prerun.in control both diagnostic output files:

```text
output_time_step_hour = 1d0
output_early_time_step_s = 300d0
output_early_until_hour = 1d0
output_wavelength_step_nm = 5d0
tau_floor = 1d-6
output_end_day = -1d0
```

This stores every 5 minutes during the first hour, then hourly. Setting
the early step to 0 disables extra samples. The main interval may instead
be specified with output_time_step_s or output_time_step_day; this also
sets the displayed time unit. Output times can be nonuniform.

For a fixed range use output_end_hour=72, or output_end_day=3.
Negative end selects automatic termination. Automatic termination waits
until all events finish and the remaining ash + formed SSA + potential
cannot exceed tau_floor at ANY model wavelength. The end is searched in
daily increments. It is a numerical negligibility threshold, not a measured
local-impact lifetime. A non-decaying case requires an explicit output end;
failure to find an end within 100 years raises an error rather than truncating
a still-active plume.

volcano_optical_depth.dat has one row for every time/layer pair and separate
total, ash and SSA cumulative optical depths for each selected wavelength.
All layers including the top contribute. Precursor is excluded.
volcano_ash_profile.dat adds an explicitly named
sulfate_precursor_potential diagnostic to its existing six particle fields.
This diagnostic is NOT optical depth or particle concentration.
No SO2 source-rate field is included in either pre-run file.

## Full Chemistry

The test holds the volcano clock until the background passes a per-species,
per-layer daily-change check for 30 consecutive days. The relative tolerance
is 1e-4 with an air-relative abundance floor of 1e-15. These are explicit
numerical criteria, not observational validation. If the default 40-year spin-up
limit is insufficient, the run stops without starting the eruption.

volcano_run.in preserves the default 730-day volcanic phase. The optional
volcano_history.in defaults to SO2 and all reactions hourly for 72 hours,
and all-wavelength solar flux every 6 hours for 72 hours. Set an output
interval <= 0 to disable that group; windows must fit inside the model run.
Run ./test only for the full spin-up plus chemistry simulation. These files
use Fortran namelists; retain the &group line and closing slash.

Formal outputs are volcano_history_species.dat, volcano_history_reactions.dat,
volcano_history_solar.dat and their reaction_map/wavelengths companions.
Solar metadata now contain a distinct integration width for EVERY bin.
The initial/final volcano_state files remain diagnostic snapshots, not restart
files. Input regeneration and make do not refresh any simulation results.
Outputs made before the 2026-09-17 spectral correction must not be treated as
results from the corrected photochemistry; a fresh spin-up and run are needed.

During volcanic evolution, the test limits outer steps, resolves event
boundaries and updates photolysis on each call. Background-only refresh
behavior is retained, except rates now initialize on the first call.
A complete long chemistry simulation and chemical time-convergence study
are still required before interpreting sulfur chemistry results.

## Viewer

The desktop program and tools/volcano_optical_depth_viewer.py read the new
and previous wide formats. Start with volcano_ash_profile.dat to inspect
local particles. The precursor field is explicitly labeled not opacity.
The default plot uses absolute values, not normalized profiles.

Column floor defaults to 1e-6 and is applied BEFORE any optional profile
normalization. Set it to 0 to inspect all nonzero values. Values below the
color minimum are white instead of being colored like significant layers.
Normalization is restricted to local fields, not cumulative shielding.
Out-of-range requests and incomplete files produce errors rather than
fabricated or silently dropped data. The stored full time range is unchanged
by display thresholds. GIF/MP4 frames represent stored samples, so extra
early frames intentionally emphasize the eruption stage.

## Regression Checks

From the repository root after make:

```bash
python3 tests/volcano_pinatubo_1991/test_prerun.py
python3 tests/volcano_pinatubo_1991/test_generation.py
```

For the independent transport checks, run make check_volcano_transport in
build. These cover conservation, positivity, bottom loss, density-weighted
equilibrium, nonuniform volumes, analytic diffusion, exponential formation,
and grid convergence. The Python checks cover formation/settling separation,
continuous source timing, full spectral/time coverage, auto-end, output
sampling, the top layer, display thresholds, malformed files and image/GIF
exports.

make check_volcano_radiation in build tests per-bin integration and J/flux
closure without running chemical time integration. test_generation.py also
checks input synchronization and the diagnostic scheduler in an isolated
fixture, not a physically equilibrated atmosphere.

check_equilibrium.f90 is invoked by test_generation.py in a temporary working
directory. It checks the consecutive-day counter, trace-abundance floor,
clock hold/reset and the failure path that must not start volcanic forcing.

For a before/after consolidation comparison, keep a separately compiled copy
of the original source tree and run:

```bash
python3 tests/volcano_pinatubo_1991/compare_consolidation.py \
  --baseline-root /path/to/original/checkout --report /tmp/volcano-comparison.json
```

The comparison runs isolated optical and frozen-atmosphere diagnostic fixtures;
it does not run the formal spin-up driver. An optional --include-chemistry adds
a short chemical integration fixture with a 180-second wall-time limit per run.
Headers, columns and output times must agree. Numerical comparisons target
byte-identical files, with a maximum relative tolerance of 5e-13 if compiler
optimization changes the last output digits. Existing model outputs are not
overwritten by these checks.
