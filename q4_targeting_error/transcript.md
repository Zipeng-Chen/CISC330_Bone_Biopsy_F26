
# Q4: Targeting error simulation - Student and Chatbot transcript

[Q4 explanation and results](README.md) | [Q1 transcript](../q1_tool_calibration/transcript.md) | [Q2 transcript](../q2_registration/transcript.md) | [Q3 transcript](../q3_biopsy_plan/transcript.md)

## 1. Perpendicular targeting error versus the 100 mm endpoint

**Student:** Does Q4 measure the distance to the needle's endpoint after 100 mm insertion, or the perpendicular distance to its trajectory?

**Chatbot:** Q4 uses the perpendicular distance from the displaced target to the registered needle axis. The insertion depth of 100 mm defines the nominal target position in the simulated scene; it is not an instruction to calculate endpoint error.

Let `P` be the registered tip, `u` its unit insertion direction, and `T'` the displaced target. Define `w=T'-P`. Then:

```text
a = dot(w, u)
closest point on needle line = P + a*u
perpendicular targeting error = norm(w - a*u)
                             = norm(cross(w, u))
```

For comparison, the endpoint after a 100 mm insertion is `P+100*u`, so:

```text
endpoint error^2 = perpendicular error^2 + (a-100)^2
```

For a hand calculation, suppose `P=(3,0,0)` mm, `u=(0,-1,0)`, and `T'=(0,-104,0)` mm. Then `w=(-3,-104,0)` and `a=104` mm. The perpendicular component is `(-3,0,0)`, so the targeting error is **3 mm**. The 100 mm endpoint is `(3,-100,0)` and is **5 mm** from the target. This distinction matters because the assignment explicitly excludes insertion-depth and core-length variation from the Q4 targeting-error metric.

The implementation evaluates this distance through `distance_point_to_line()` in [simulation.py](simulation.py). It treats the trajectory as a line; whether the closest point lies along the physically reachable forward ray is a separate question.

## 2. Constructing the ideal scene and checking the noise-free result

**Student:** How do we choose the Q4 geometry, and what should happen before adding noise?

**Chatbot:** The model and CT coordinate frames coincide for Q4. The implementation uses millimetres and the assignment's tool-marker geometry:

| Quantity | Coordinates |
| --- | --- |
| Marker 1 | `(20,180,0)` mm |
| Marker 2 | `(80,180,0)` mm |
| Marker 3 | `(20,240,0)` mm |
| Tool tip | `(0,0,0)` mm |
| Needle axis | `(0,-1,0)` |
| Planned target | `(0,-100,0)` mm |

The target is exactly 100 mm ahead of the ideal tip. With no localization noise or displacement, registration should recover the original trajectory and the perpendicular targeting error should be zero, apart from floating-point roundoff.

A direct hand check gives `w=(0,-100,0)` and `dot(w,u)=100`. Therefore:

```text
w - dot(w,u)*u
= (0,-100,0) - 100*(0,-1,0)
= (0,0,0)

error = 0 mm
```

The numerical noise-free check gave approximately `9.89e-14` mm, which is effectively zero. This is a useful sanity check before interpreting Monte Carlo results. It does not by itself establish that the noisy experiment is correct.

## 3. Converting the tracker's 95% radial bound into Gaussian noise

**Student:** The assignment says the tracker error magnitude is less than 1 mm in 95% of observations. Why can't we just use `sigma=1` mm for each coordinate?

**Chatbot:** The statement constrains the **3D error magnitude**, not the standard deviation of each coordinate. If the tracker noise is isotropic Gaussian, write:

```text
ex, ey, ez independently follow Normal(0, sigma_tracker^2)
radial error R = sqrt(ex^2 + ey^2 + ez^2)
P(R <= 1 mm) = 0.95
```

Here `R/sigma_tracker` follows a chi distribution with three degrees of freedom. Solving its 95th-percentile condition gives:

```text
sigma_tracker ≈ 0.357719874 mm per coordinate
```

Using 1 mm as the coordinate sigma would produce a much larger radial error and would not represent the stated tracker accuracy. An independent 100,000-sample check produced about **95.02%** of radial errors within 1 mm, consistent with the requested 95% bound.

The code reuses the shared `sigma_for_radial_bound()` helper. Each CT marker receives an independent 3D Gaussian perturbation with this fixed coordinate sigma. The tracker sigma remains constant during the procedural-displacement sweep.

## 4. Propagating marker uncertainty through registration

**Student:** If the target displacement is zero, why can the needle still miss by more than 3 mm?

**Chatbot:** Zero target displacement does **not** mean zero registration error. The three marker locations observed in CT are noisy. Registering the model markers to those observations estimates a rigid transformation whose translation and rotation can differ from the ideal transformation. The resulting tip and axis are therefore displaced or tilted.

In the primary experiment, a trial performs the following steps:

1. Start with the exact tool tip, axis, markers, and planned target.
2. Perturb the three CT marker positions independently using fixed tracker Gaussian noise.
3. Estimate the tool-to-CT registration from the model markers and noisy CT markers.
4. Draw an independent 3D target displacement with coordinate standard deviation `sigma_procedural`.
5. Measure the perpendicular distance between the displaced target and registered needle line.
6. Count the trial as successful if that distance is at most 3 mm.

For intuition, even a small orientation error can produce substantial lateral error 100 mm away from the tip. A direction tilted by `theta` radians contributes approximately `100*sin(theta)` mm of lateral displacement at the planned depth, before considering tip translation. This is an illustration of error amplification, not a replacement for the full registration simulation.

A simplified hand example is a line with tip `(0,0,0)` and unit direction `(sin(theta),-cos(theta),0)`, aimed at the fixed target `(0,-100,0)`. Its perpendicular target distance is `100*abs(sin(theta))`. At `theta=2°`, that is about **3.49 mm**, already above the tolerance. In the actual experiment, both tip and direction can change, so the code evaluates the complete point-to-line distance on every trial.
