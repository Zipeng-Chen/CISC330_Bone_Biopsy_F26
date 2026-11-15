# Q2: Registration - Student and Chatbot transcript


[Q2 explanation and results](README.md) | [Q1 transcript](../q1_tool_calibration/transcript.md) | [Q3 transcript](../q3_biopsy_plan/transcript.md)

## 1. Figure 2 units and the shared model

**Student:** Why are Figure 2 coordinates multiplied by 10, and which quantities must remain unchanged?

**Chatbot:** The professor confirmed that Figure 2's coordinates are in cm.
The code uses mm, so it converts those positions by multiplying by 10. This
is a unit conversion, not a change to the physical tool size.

| Marker | Figure 2 coordinates (cm) | Coordinates used by the code (mm) |
| --- | --- | --- |
| M1 | (2,18,0) | (20,180,0) |
| M2 | (8,18,0) | (80,180,0) |
| M3 | (2,24,0) | (20,240,0) |

The triangle has two 60 mm legs. Treating the original numbers as mm would
instead give 6 mm legs and a much smaller tool. With the same absolute marker
noise, that smaller layout gives much greater orientation uncertainty.

`assignment_tool_model()` in [utils/models.py](../utils/models.py) already
performs this conversion. Q1, Q2, Q3, and Q4 reuse that model, so its returned
coordinates must not be multiplied by 10 again.

Any nonzero tip offset or other position originally expressed in cm must
also be converted to mm. An offset already supplied in mm stays unchanged.
The assignment's tip is `(0,0,0)`, so converting it leaves it at zero.

These quantities do not change:

- The insertion unit direction `(0,-1,0)`, which has no length unit.
- Rotation matrices and rotation angles.
- The marker triangle's shape and relative proportions.
- Values already specified in mm: the 1 mm marker-error bound, 3 mm targeting
  tolerance, 100 mm target depth, and any CT tumor center or radius already
  expressed in mm.

All positional inputs to registration must use the same unit. The algorithm
fits only rotation and translation; it cannot repair mismatched units by
resizing one set of markers.

## 2. Recovering the least-squares rotation with SVD

**Student:** How does SVD recover the least-squares rotation?

**Chatbot:** Given corresponding model points `a_i` and CT points `b_i`, the
rigid-registration problem is to minimize:

```text
E = sum_i norm(R*a_i + t - b_i)^2
```

over a proper rotation `R` and translation `t`. Let `a_bar` and `b_bar` be
the two centroids. Minimizing with respect to translation gives:

```text
t = b_bar - R*a_bar.
```

After subtracting the centroids, write `a_i' = a_i-a_bar` and
`b_i' = b_i-b_bar`. Substituting and expanding gives:

```text
E = sum_i norm(a_i')^2 + sum_i norm(b_i')^2
    - 2*sum_i (b_i')^T R a_i'.
```

The first two terms do not depend on `R`. Minimizing the error therefore
means maximizing the final inner-product sum. Define:

```text
H = sum_i a_i' (b_i')^T
  = X.T @ Y
```

where `X` and `Y` store centered points as rows. The sum to maximize is
`trace(R H)`. With SVD `H = U S V.T`:

```text
trace(R H) = trace(S V.T R U).
```

Set `M = V.T R U`. Since `M` is orthogonal, each diagonal entry is at most
one. With nonnegative singular values `s_i`:

```text
trace(S M) = sum_i s_i*M_ii <= sum_i s_i.
```

For an unconstrained orthogonal transformation, the bound can be reached by
`M=I`, giving `R=V U.T`. However, an orthogonal matrix can be a reflection.
A physical rotation must have determinant `+1`.

If `det(V U.T)=-1`, flip the direction associated with the smallest singular
value. This gives the smallest loss in the objective while enforcing a
proper rotation:

```text
D = diag(1, 1, det(V U.T))
R = V D U.T
t = b_bar - R*a_bar.
```

Three noncollinear markers are sufficient even though they are coplanar.
Their centered configuration has rank two, so a zero third singular value
is expected. The proper-rotation constraint fixes the handedness of the
third direction. Collinear markers are different: rotation about their line
remains unresolved, so the code rejects that case.

The implementation is `register_rigid()` in
[registration.py](registration.py). It checks the inputs, fits this proper
rotation, recovers translation, and reports the RMS marker-fit residual.
The three exact translation/rotation cases and their hand predictions are
documented in the [Q2 README](README.md).

## 3. Transforming points and directions

**Student:** Why does the tip receive rotation and translation, while the axis receives rotation only?

**Chatbot:** A point has a location, so translating the tool moves that point.
A direction describes orientation, so moving the whole tool without turning
it cannot change the direction.

Take two points on the needle, `P` and `Q`, with direction proportional to
their difference `v=Q-P`. Under the rigid transformation:

```text
P_CT = R P + t
Q_CT = R Q + t
Q_CT - P_CT = R(Q-P) + (t-t) = R v.
```

Translation cancels. Therefore:

```text
tip_CT  = R @ tip_model + t
axis_CT = R @ axis_model.
```

In homogeneous coordinates, a point has fourth component 1 and a direction
has fourth component 0. Multiplying by `[R t; 0 1]` adds translation only to
the point.

The required test cases make this distinction visible:

| Transformation | Tip result | Axis result |
| --- | --- | --- |
| Translation by (10,-5,2) mm | Moves to (10,-5,2) mm | Remains (0,-1,0) |
| +90 degrees about z | Origin remains the origin | Becomes (1,0,0) |
| That rotation, then the translation | Moves to (10,-5,2) mm | Remains (1,0,0) |

A proper rotation preserves unit length. The implementation still normalizes
direction inputs and outputs so a valid nonunit input has the same meaning.
Checking only the output norm is not enough to detect a bad transform: an
incorrectly translated direction could be normalized and still point the
wrong way. The tests compare the expected direction itself and verify that
pure translation leaves it unchanged.

The zero tip does not move under rotation about the origin. For that reason,
the tests also use a nonzero tip and an out-of-plane rotation, so an incorrect
tip transformation cannot be hidden by that special case.

## 4. Marker-fit residual versus accuracy at the tip

**Student:** Why does a small marker-fit residual not necessarily mean a small error at the needle tip?

**Chatbot:** The residual measures agreement at the markers after the fitted
transformation. The tip error measures agreement at a different location,
using the true tip as the reference. They answer different questions.

For the RMS marker residual:

```text
residual_i = R @ model_marker_i + t - observed_marker_i
RMS = sqrt(mean(norm(residual_i)^2)).
```

The fit has already chosen `R` and `t` to reduce these residuals. Some
measurement error is absorbed into the estimated pose instead of remaining
visible in the residual. A low residual therefore does not imply that the
estimated pose equals the true pose.

A simple separate illustration is a common shift of all observed markers
by `(3,0,0)` mm. Registration can fit that shift exactly, leaving zero marker
residual, while the estimated tip is displaced by 3 mm from the true tip.
This example is a coherent measurement offset; it is not the independent
Gaussian error model used in the assignment simulation.

Rotation errors also grow in their effect away from the marker cluster.
For a small angular error `delta_theta` in radians, the displacement due to
rotation is approximately:

```text
error_from_rotation = delta_theta cross r
```

where `r` is measured from the rotation reference point. The magnitude depends
on the perpendicular distance to that rotation axis. A 0.5-degree angular
error over an effective 300 mm distance gives about
`300 * (0.5*pi/180) = 2.62` mm of displacement.

For this tool, the marker centroid is `(40,200,0)` mm. The tip is at the
origin, and the ideal 100 mm target is at `(0,-100,0)`. Both lie outside the
marker triangle. The triangle's spread controls how well orientation can be
estimated; distance from that cluster then affects errors at the tip and
target. Translation and rotation errors from the same fit can also be
correlated, so their magnitudes should not simply be added.

The marker-fit residual is a useful diagnostic, often called fiducial
registration error. In simulation, the useful accuracy checks also include
tip error and axis-angle error against ground truth, and the actual
point-to-line targeting error. A small residual alone is not a targeting
guarantee.

The noise should enter at the marker observations before this fit. That
procedure is explained in the [Q1 transcript](../q1_tool_calibration/transcript.md).
The distinction between needle-axis and endpoint distances is explained in
the [Q3 transcript](../q3_biopsy_plan/transcript.md).
