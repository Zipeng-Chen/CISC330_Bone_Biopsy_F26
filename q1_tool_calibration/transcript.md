# Q1: Tool calibration - Student and Chatbot transcript

[Q1 explanation and results](README.md) | [Q2 transcript](../q2_registration/transcript.md) | [Q3 transcript](../q3_biopsy_plan/transcript.md)

## 1. Pivot motion and independent orientations

**Student:** Why must pivoting include sufficiently different orientations about multiple axes?

**Chatbot:** Pivot calibration uses the fact that the physical tip stays at
one fixed point while the marker frame moves around it. The unknowns are the
tip position `p` in the model frame and the fixed pivot `q` in the tracker
frame. Each registered pose must satisfy:

```text
R_i p + t_i = q
[R_i  -I] [p; q] = -t_i.
```

Stacking all poses gives a system with six unknown coordinates. The solution
is unique only when the stacked matrix has rank six. More equations alone
are not enough if they repeat the same geometric information.

To see the problem with a single rotation axis, first consider aligned
frames and an axis `a` that every rotation leaves unchanged: `R_i a = a`.
If `(p,q)` is a solution, then:

```text
R_i (p + lambda*a) + t_i = q + lambda*a
```

also holds for any `lambda`. Moving the estimated tip along that axis cannot
be distinguished from moving the estimated pivot along it. The system has
rank at most five; repeated identical orientations leave even more freedom.
With a different initial tracker orientation, the invariant model direction
maps to one constant tracker direction, and the same ambiguity remains.

Physically, spinning a shaft cannot identify where along the shaft its tip
is. Sufficiently varied tilts about different axes remove this ambiguity.
Two poses alone are still insufficient: their relative rotation leaves an
axis unresolved. The implementation requires at least three poses and checks
the actual rank rather than relying on a pose count alone.

With noisy observations, conditioning matters as well. Very small tilts can
produce a full-rank system that is still sensitive to measurement error.
The poses should cover different tilt directions with enough angular spread
to provide useful information. The current simulation uses 60 pivot poses
within a cone with a 40-degree half-angle. These are chosen experiment
settings, not an assertion that this motion is optimal.

The simulated translation is `t_i = q - R_i @ true_tip`, which keeps the tip
fixed. The calibration solver receives the marker observations and recovers
the poses itself; it does not receive the true tip or true poses.

The code reports a condition number and rejects rank-deficient systems.
Those checks help identify poor motion, but a small residual alone does not
prove an accurate tip. See [calibrate_tip](calibration.py) and the motion
checks in [test_calibration.py](test_calibration.py).

## 2. Why a fixed-sleeve spin produces circles

**Student:** Why does spinning the shaft in a stationary sleeve make off-axis markers trace planar circles?

**Chatbot:** The sleeve holds one axis line fixed while the shaft rotates
about it. Let `c` be a point on that line and `a` its unit direction. Split a
marker's initial position into an axial part and a perpendicular part:

```text
m_0 = c + h*a + rho
h = dot(m_0 - c, a)
dot(rho, a) = 0.
```

A rotation about `a` leaves the axial part unchanged and rotates the radial
part:

```text
m(theta) = c + h*a + R(theta)*rho.
```

Two quantities stay constant. First, `dot(m(theta)-c, a)=h`, so the marker
stays in a plane perpendicular to the shaft. Second, rotation preserves
`norm(rho)`, so the distance from the circle center `c+h*a` stays fixed.
Together, these constraints describe a circle.

In the model frame, the shaft is on the y axis. With the confirmed mm model:

| Marker | Circle center (mm) | Radius (mm) | Plane |
| --- | --- | ---: | --- |
| M1 | (0,180,0) | 20 | y = 180 |
| M2 | (0,180,0) | 80 | y = 180 |
| M3 | (0,240,0) | 20 | y = 240 |

All three planes have normals parallel to the shaft. M1 and M2 have the same
circle center, but different radii. A marker exactly on the shaft would have
zero radius and would not supply a plane normal.

The implementation fits one plane per usable marker with SVD, aligns the
normal signs, and combines the normals using weights based on their in-plane
spread. It rotates the tracker-frame direction back into the model frame
using `R_i.T`. It does not translate a direction.

Circle fitting and a line through distinct fitted circle centers would be
another possible approach. Our implementation uses plane normals; it does
not claim to perform that extra circle-center fit.

The movement must avoid both tilting and sliding along the sleeve. Sliding
can turn the circle into a helix, while tilting changes the plane. A short
arc gives little spread relative to the noise; a full turn provides broader
coverage. The current simulation uses 60 positions across a complete turn.
The data-quality checks reject unresolved or inconsistent planes, but they
cannot guarantee that every physical constraint was followed perfectly.

## 3. Selecting the insertion-facing sign

**Student:** How does the calibrated tip resolve the insertion-facing sign, and what assumption about marker placement is required?

**Chatbot:** A plane fit gives a normal only up to sign. Both `a` and `-a`
describe the same plane, so spinning alone does not identify the insertion
direction.

After bringing the axis into the model frame, form the vector from the model
marker centroid to the calibrated tip:

```text
toward_tip = calibrated_tip - mean(model_markers).
```

If `dot(a, toward_tip)` is negative, reverse `a`. If it is positive, keep it.
The chosen unit direction then points from the marker side toward the tip.
If the dot product is too close to zero, the sign is unresolved and the code
raises an error.

For this model, the marker centroid is `(40,200,0)` mm and the ideal tip is
`(0,0,0)`. Therefore `toward_tip=(-40,-200,0)`. A fitted normal `(0,1,0)` has
dot product `-200` and is reversed to `(0,-1,0)`.

The physical assumption is that the marker cluster is behind the working
tip along the insertion direction. If markers were placed beyond the tip,
the same rule would choose the wrong sign. If their centroid were at the
same axial level as the tip, the dot product could be nearly zero and noise
could make its sign unstable. The radial marker offset does not decide this
sign; only the projection onto the shaft direction matters.

The rule uses the estimated tip from pivot calibration, not the known true
axis. The large axial separation makes it relatively insensitive to small
tip errors for this particular tool. It is still a tool-design assumption,
not a universal method for every marker layout.

## 4. Understanding the 1 mm probability requirement

**Student:** Why does “95% of marker errors are below 1 mm” refer to the 3D error magnitude, rather than each coordinate separately?

**Chatbot:** The requirement concerns the displacement of a marker in space.
Its magnitude is:

```text
norm(e) = sqrt(e_x^2 + e_y^2 + e_z^2).
```

This length is unchanged if the coordinate axes are rotated. A bound on one
coordinate describes a different quantity.

Let the three coordinate errors be independent Gaussian variables with mean
zero and standard deviation `sigma`. For a single coordinate, the familiar
95% interval would give `sigma = 1/1.96`, approximately 0.5102 mm. However,
that would put only about 72% of the three-dimensional error vectors inside
a 1 mm sphere. About 28%, rather than 5%, would lie outside it.

For the vector magnitude, `norm(e)/sigma` has a chi distribution with three
degrees of freedom. This gives approximately **0.35772 mm per coordinate**
for the required 95% radial bound. The one-dimensional `1/1.96` choice does
not meet that requirement for the full displacement vector.

The simulation draws independent errors per coordinate, marker, and pose.
It does not clip the magnitude at 1 mm, since larger errors must remain
possible under the Gaussian model.

## 5. Deriving the coordinate standard deviation

**Student:** How does that requirement produce approximately 0.35772 mm standard deviation per coordinate?

**Chatbot:** Divide the sum of the three squared coordinate errors by their
common variance:

```text
X = (e_x^2 + e_y^2 + e_z^2) / sigma^2
X follows chi-square with 3 degrees of freedom.
```

The requirement is `P(norm(e) < 1 mm) = 0.95`. Therefore:

```text
P(X < (1 mm / sigma)^2) = 0.95
(1 mm / sigma)^2 = 7.814727903
sigma = 1 mm / sqrt(7.814727903)
      = 0.357719874 mm per coordinate.
```

The three different error summaries should not be mixed up:

| Quantity | Value |
| --- | ---: |
| Coordinate standard deviation | 0.35772 mm |
| 3D RMS marker error, sqrt(3)*sigma | About 0.620 mm |
| 95th percentile of 3D marker-error magnitude | 1.000 mm |

The project obtains `sigma` with `sigma_for_radial_bound()` in
[utils/noise.py](../utils/noise.py). It uses bisection on the radial CDF:

```text
F(s) = erf(s/sqrt(2)) - sqrt(2/pi)*s*exp(-s^2/2)
F(1 mm / sigma) = 0.95.
```

This computes the same value without adding SciPy as a dependency. A large
random sample should have approximately 95% of its magnitudes below 1 mm;
the exact empirical percentage varies between samples.

The saved Q1 experiment contains 180,000 marker observations from 500 noisy
calibrations, with 60 poses for each of the two movements. Its measured
fraction below 1 mm is **95.0439%**, consistent with the specified model.
The settings and results are in
[calibration-results.json](outputs/calibration-results.json).

The target-displacement sigma swept in Q4 is also a per-coordinate standard
deviation, but it is a separate parameter. It is not substituted for the
fixed instrument-noise sigma.

## 6. Adding noise before calibration and registration

**Student:** Why must noise be added to marker observations before registration, rather than directly to the final calibrated tip?

**Chatbot:** The measurement error enters at the markers. The fitted rotation
and translation, and then the estimated tip and axis, are calculated from
those observations. Adding arbitrary noise to a final tip bypasses the
process whose accuracy the experiment is meant to measure.

For Q1, the sequence is:

1. Generate exact marker observations for the constrained pivot and spin
   movements.
2. Add independent Gaussian errors to the observed marker coordinates.
3. Recover each pose with `register_rigid()` using the noisy observations.
4. Solve the pivot and plane-normal calibration problems again.
5. Compare the resulting tip and axis with the simulated ground truth.

The marker layout matters because a small rotation error can move a distant
tip substantially. The fitted pose errors can be correlated, and the
resulting tip and axis errors need not be isotropic or independent. A single
hand-chosen tip-noise distribution cannot capture those effects without a
separate justification.

The same principle applies to CT registration in Q4: add noise to the CT
markers, re-register, and transform the calibrated tip and axis using that
same fitted pose. If calibration uncertainty is included too, keep each
trial's estimated tip and axis together and use independent CT measurement
noise. Independently perturbing an already computed tip and axis would lose
their relationship.

Procedural target displacement is different: it represents movement of the
target, so it is added to the target point independently of the instrument
observations. It is not added to the markers.

The Q1 integration example uses exact calibrated values to check the known
geometry. The separate noisy experiment measures calibration accuracy. Those
two checks should not be described as a full noisy targeting simulation.

