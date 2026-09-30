"""Constant-D diffusion into a slab exposed to equal surface reservoirs.

The solver uses z=x/L, tau=D*t/L**2, and u=C/C_surface. Both end
points are fixed at u=1 from tau=0 onward; initially the interior is empty.
See docs/METHOD.md for the model, assumptions, and series derivation.
"""

from dataclasses import dataclass
from numbers import Integral

import numpy as np


@dataclass(frozen=True)
class Snapshot:
    tau: float
    position: np.ndarray
    concentration: np.ndarray

    @property
    def uptake(self) -> float:
        """Trapezoidal estimate of the normalized slab-average concentration."""
        return float(np.trapezoid(self.concentration, self.position))


def solve_slab(times, nodes: int = 101, r: float = 0.45) -> list[Snapshot]:
    """Solve u_tau=u_zz with FTCS; return independent snapshots at `times`.

    `times` are strictly increasing, finite, non-negative Fourier times.
    `r` is the largest permitted dtau/dz**2 and must lie in (0, 0.5].
    Steps shorten to land exactly on each requested output time.
    At tau=0, end points are already set to one (the boundary jump at 0+).
    """
    if isinstance(nodes, bool) or not isinstance(nodes, Integral) or nodes < 3:
        raise ValueError('nodes must be an integer >= 3')
    if not np.isfinite(r) or not 0 < r <= 0.5:
        raise ValueError('FTCS requires 0 < r <= 0.5')
    targets = np.asarray(times, dtype=float)
    if (targets.ndim != 1 or targets.size == 0
            or not np.all(np.isfinite(targets)) or np.any(targets < 0)
            or np.any(np.diff(targets) <= 0)):
        raise ValueError('times must be finite, non-negative, and strictly increasing')

    position = np.linspace(0, 1, nodes)
    dz = 1 / (nodes - 1)
    base_step = r * dz**2
    concentration = np.zeros(nodes)
    concentration[[0, -1]] = 1
    tau = 0.0
    output = []
    for target in targets:
        remaining = float(target) - tau
        # Split the interval into identical stable steps; this avoids a
        # nearly zero final step caused by floating-point accumulation.
        steps = int(np.ceil(remaining / base_step))
        if steps:
            ratio = (remaining / steps) / dz**2
            for _ in range(steps):
                concentration[1:-1] += ratio * (
                    concentration[:-2] - 2 * concentration[1:-1]
                    + concentration[2:]
                )
        tau = float(target)
        output.append(Snapshot(tau, position.copy(), concentration.copy()))
    return output


def _modes(tau: float, terms: int) -> np.ndarray:
    if not np.isfinite(tau) or tau < 0:
        raise ValueError('tau must be finite and non-negative')
    if isinstance(terms, bool) or not isinstance(terms, Integral) or terms < 1:
        raise ValueError('terms must be a positive integer')
    return np.arange(1, 2 * terms, 2, dtype=float)


def analytical_profile(position, tau: float, terms: int = 512) -> np.ndarray:
    """Odd sine series for this slab problem, with a finite term count.

    The study checks times >= 0.002; extremely small positive times may
    need more terms to resolve the initial boundary discontinuity.
    """
    odd = _modes(tau, terms)
    z = np.asarray(position, dtype=float)
    if z.ndim != 1 or not np.all(np.isfinite(z)) or np.any((z < 0) | (z > 1)):
        raise ValueError('position must be a finite 1D array in [0, 1]')
    boundary = (z == 0) | (z == 1)
    if tau == 0:
        return boundary.astype(float)
    amplitudes = np.exp(-(odd * np.pi)**2 * tau) / odd
    result = 1 - (4 / np.pi) * (amplitudes @ np.sin(np.pi * odd[:, None] * z))
    result[boundary] = 1
    return result


def analytical_uptake(tau: float, terms: int = 512) -> float:
    """Exact spatial integral of the truncated concentration series."""
    odd = _modes(tau, terms)
    if tau == 0:
        return 0.0
    return float(1 - 8 / np.pi**2 * np.sum(np.exp(-(odd*np.pi)**2*tau) / odd**2))


def time_for_uptake(fraction: float = 0.9) -> float:
    """Fourier time for 1% <= uptake < 100%, by bisection of the series.

    Smaller fractions are excluded because this finite series is not
    intended to resolve the extremely early-time boundary layer.
    """
    if not np.isfinite(fraction) or not 0.01 <= fraction < 1:
        raise ValueError('fraction must satisfy 0.01 <= fraction < 1')
    low, high = 0.0, 1.0
    while analytical_uptake(high) < fraction:
        high *= 2
    for _ in range(80):
        mid = (low + high) / 2
        if analytical_uptake(mid) < fraction:
            low = mid
        else:
            high = mid
    return (low + high) / 2


def physical_time(tau, thickness_m: float, diffusivity_m2_s: float):
    """Convert Fourier time to seconds using t=tau*L**2/D (SI units)."""
    if (not np.isfinite(thickness_m) or thickness_m <= 0
            or not np.isfinite(diffusivity_m2_s) or diffusivity_m2_s <= 0):
        raise ValueError('thickness and diffusivity must be finite and positive')
    values = np.asarray(tau, dtype=float)
    if not np.all(np.isfinite(values)) or np.any(values < 0):
        raise ValueError('tau must be finite and non-negative')
    return values * thickness_m**2 / diffusivity_m2_s
