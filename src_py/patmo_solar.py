"""Sample photon spectral flux on the actual generated PATMO bin edges."""
from bisect import bisect_left, bisect_right
from math import isfinite
from pathlib import Path
from patmo_numeric import linear_interpolator

HC_EV_NM = 4.135667662e-15 * 2.99792458e10 * 1e7


def bin_average(wavelength, flux, energy_left, energy_right):
    wavelength, flux = list(map(float, wavelength)), list(map(float, flux))
    left, right = list(map(float, energy_left)), list(map(float, energy_right))
    if len(wavelength) < 2 or len(flux) != len(wavelength):
        raise ValueError("Solar spectrum requires at least two wavelength/flux pairs")
    if not all(map(isfinite, wavelength)) or any(b <= a for a, b in zip(wavelength, wavelength[1:])):
        raise ValueError("Solar wavelengths must be finite, increasing and unique")
    if any(not isfinite(f) or f < 0 for f in flux):
        raise ValueError("Photon flux must be finite and nonnegative")
    if not left or len(left) != len(right):
        raise ValueError("Invalid photochemistry bin arrays")
    if any(not isfinite(a+b) or a <= 0 or b <= a for a, b in zip(left, right)):
        raise ValueError("Invalid photochemistry energy edges")
    low, high = [HC_EV_NM/e for e in right], [HC_EV_NM/e for e in left]
    tolerance = 1e-10 * max(max(high), 1)
    if wavelength[0] > min(low)+tolerance or wavelength[-1] < max(high)-tolerance:
        raise ValueError("Solar spectrum does not cover all model wavelength bins")
    interpolate = linear_interpolator(wavelength, flux)
    result = []
    for lo, hi in zip(low, high):
        x = [lo] + wavelength[bisect_right(wavelength, lo):bisect_left(wavelength, hi)] + [hi]
        y = [interpolate(value) for value in x]
        integral = sum(0.5*(a+b)*(v-u) for a, b, u, v in zip(y, y[1:], x, x[1:]))
        result.append(integral/(hi-lo))
    return result


def numeric_rows(path):
    with open(path) as stream:
        return [[float(value) for value in line.split("#", 1)[0].split()]
                for line in stream if line.split("#", 1)[0].strip()]


def build_solar_flux(test_name, metric_path="build/xsecs/photoMetric.dat"):
    source = Path("tests") / test_name
    workbook = source / "solar_flux.xlsx"
    metric = numeric_rows(metric_path)
    if not workbook.exists():
        values = numeric_rows("build/solar_flux.txt")
        if len(values) != len(metric) or any(len(row) != 1 or not isfinite(row[0]) or row[0] < 0 for row in values):
            raise ValueError("Invalid text-only solar flux: require one nonnegative value per bin")
        print("Using supplied solar_flux.txt; caller is responsible for bin order and photon/nm units")
        return
    try:
        from openpyxl import load_workbook
    except ModuleNotFoundError as exc:
        if exc.name != "openpyxl":
            raise
        raise RuntimeError("Reading solar_flux.xlsx requires openpyxl. Cannot safely skip solar-spectrum conversion.") from None
    book = load_workbook(workbook, read_only=True, data_only=True)
    try:
        rows = iter(book.worksheets[0].values)
        header = next(rows)
        wcol = next(i for i, c in enumerate(header) if str(c).strip().lower().startswith("wavelength"))
        fcol = next(i for i, c in enumerate(header) if str(c).strip().lower().startswith("irradiance"))
        pairs = sorted((float(row[wcol]), float(row[fcol])) for row in rows
                       if row[wcol] is not None or row[fcol] is not None)
    finally:
        book.close()
    x, y = zip(*pairs)
    values = bin_average(x, y, [row[2] for row in metric], [row[3] for row in metric])
    for path in (source / "solar_flux.txt", Path("build/solar_flux.txt")):
        path.write_text("".join("%.15e\n" % value for value in values))
    print("Solar photon spectrum aligned to", len(values), "actual PATMO bins (long to short wavelength)")
