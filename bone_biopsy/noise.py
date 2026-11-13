"""Isotropic Gaussian marker noise with a specified 3D radial probability."""

from math import erf, exp, pi, sqrt

import numpy as np

from ._validation import nonnegative_scalar


def sigma_for_radial_bound(bound_mm=1.0, probability=0.95):
    """Find coordinate sigma such that P(||N(0, sigma^2 I3)|| < bound)=probability.

    The standardized 3D length follows a chi distribution with three degrees
    of freedom. A 1 mm radial 95% bound gives sigma about 0.35772 mm, not 0.5 mm.
    """
    bound = nonnegative_scalar(bound_mm, "bound_mm")
    chance = nonnegative_scalar(probability, "probability")
    if bound == 0 or not 0 < chance < 1:
        raise ValueError("bound_mm must be positive and probability must lie in (0, 1)")

    def radial_cdf(distance):
        return erf(distance / sqrt(2)) - sqrt(2 / pi) * distance * exp(-distance**2 / 2)

    # Invert the radial CDF without introducing a SciPy dependency
    lower, upper = 0.0, 1.0
    while radial_cdf(upper) < chance:
        upper *= 2
    for _ in range(80):
        middle = (lower + upper) / 2
        if radial_cdf(middle) < chance:
            lower = middle
        else:
            upper = middle
    return bound / ((lower + upper) / 2)


def add_marker_noise(observations, sigma_mm, *, rng):
    """Return new observations with independent noise per frame, marker, and axis"""
    points = np.asarray(observations, dtype=float)
    if points.ndim < 2 or points.shape[-1] != 3 or not np.all(np.isfinite(points)):
        raise ValueError("observations must be finite with last dimension 3")
    sigma = nonnegative_scalar(sigma_mm, "sigma_mm")
    if not isinstance(rng, np.random.Generator):
        raise TypeError("rng must be a numpy.random.Generator")
    # Add noise to marker measurements, before estimating any tool poses.
    return points + rng.normal(0.0, sigma, size=points.shape)
