"""Linear interpolation without third-party dependencies."""
from bisect import bisect_right
from math import isfinite


def linear_interpolator(x, y, outside=None):
    x, y = list(map(float, x)), list(map(float, y))
    if not x or len(x) != len(y):
        raise ValueError("Interpolation requires matching nonempty arrays")
    x, y = map(list, zip(*sorted(zip(x, y), key=lambda pair: pair[0])))
    if not all(map(isfinite, x + y)):
        raise ValueError("Interpolation data must be finite")

    def interpolate(value):
        if value < x[0]:
            return y[0] if outside is None else outside
        if value > x[-1]:
            return y[-1] if outside is None else outside
        i = bisect_right(x, value) - 1
        if i == len(x) - 1 or value == x[i]:
            return y[i]
        return y[i] + (value-x[i]) * ((y[i+1]-y[i]) / (x[i+1]-x[i]))
    return interpolate
