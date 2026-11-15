# Q3: Biopsy planning - Student and Chatbot transcript

[Q3 explanation and results](README.md) | [Q1 transcript](../q1_tool_calibration/transcript.md) | [Q2 transcript](../q2_registration/transcript.md)

## 1. Needle-axis distance and the 100 mm endpoint

**Student:** How does Figure 6 define targeting error, and why is it different from distance to the endpoint after 100 mm insertion?

**Chatbot:** Figure 6 shows the perpendicular distance from the needle axis
to the displaced target. The red segment meets the trajectory at a right
angle. This agrees with the professor's direction to use Figure 6 when
interpreting the targeting-error requirement.

Let `P` be the registered tip, `u` the unit insertion direction, and `T` the
displaced target. Set `w=T-P`. Then:

```text
a = dot(w,u)
closest point on the line = P + a*u
targeting error d = norm(w - a*u)
                  = norm(cross(w,u)).
```

By comparison, the endpoint after exactly 100 mm is `P+100*u`. Its distance
to the target includes both a perpendicular component and an along-axis
component:

```text
endpoint distance^2 = d^2 + (a-100)^2.
```

The endpoint distance is therefore at least as large as the point-to-line
distance. They are equal only when the target's closest point lies at the
100 mm endpoint.

For example, let the actual tip be `(3,0,0)` mm and the direction `(0,-1,0)`.
The planned target was `(0,-100,0)` mm, but it moves to `(0,-104,0)` mm.
The needle line is `(3,-s,0)`. Its closest point to the displaced target is
`(3,-104,0)`, giving **3 mm targeting error**. The 100 mm endpoint is instead
`(3,-100,0)`, whose distance to the displaced target is
`sqrt(3^2+4^2)=5` mm.

The 100 mm value sets the planned target depth in the ideal scene. It does
not change the Figure 6 metric into endpoint distance. The assignment also
sets aside uncertainty in insertion depth and core length in its Q4 error
simulation. Q3 still calculates those planning quantities separately.

`distance_point_to_line()` measures an infinite line, so it does not establish
that the target lies ahead of the current tip. Forward insertion requires
the separate checks described below.

Although the perpendicular error has two components, its magnitude is not
automatically Rayleigh-distributed. That distribution applies when those
components are independent, zero-mean Gaussian variables with equal variance.
Registration can create unequal and correlated errors even when the marker
noise is isotropic. Use the empirical 95th percentile for the full simulation
rather than assuming a fixed Rayleigh or `1.96*sigma` multiplier.

## 2. Normalizing the insertion direction

**Student:** Why must the needle direction be normalized before calculating insertion distances?

**Chatbot:** The line is written `X(s)=P+s*u`. When `u` has unit length, moving
the parameter by `s` moves the tip by `s` mm. With a nonunit vector `v`, the
physical distance is instead `abs(s)*norm(v)`.

For example, using `(0,-10,0)` directly would make `s=100` correspond to
1,000 mm of movement. The vector points in the correct direction, but its
length makes `s` an arbitrary parameter instead of an insertion distance.

The projection formula also depends on normalization. For unit `u`, the
along-axis distance is `dot(T-P,u)`. For a general nonzero `v`, the closest
point parameter would be `dot(T-P,v)/dot(v,v)`. Similarly:

```text
Unit direction:     line distance = norm(cross(T-P,u))
Nonunit direction:  line distance = norm(cross(T-P,v)) / norm(v).
```

The same issue appears in the sphere intersection. For a general direction,
substituting `P+s*v` into `norm(X-C)^2=r^2` gives:

```text
A = dot(v,v)
B = 2*dot(v,P-C)
C0 = dot(P-C,P-C) - r^2
A*s^2 + B*s + C0 = 0.
```

Those roots are parameters along `v`; their difference must be multiplied by
`norm(v)` to obtain a physical length. Normalizing first makes the returned
roots signed distances in mm and simplifies the remaining calculations.

The project normalizes directions and rejects zero or invalid vectors.
A useful check is that a deliberately nonunit input, such as `(0,0,10)`,
produces the same plan as `(0,0,1)`. This case is covered by the planning tests.

The two depths must be distinguished: the near root is the entry depth,
and the far root is the depth needed to traverse the full available core.

## 3. Inserting to the far tumor boundary

**Student:** Why does full-core insertion depth reach the far tumor boundary?

**Chatbot:** The assignment treats the available core as the forward segment
of the needle trajectory inside the spherical tumor. Reaching the near
boundary only reaches the start of that segment. Traversing the full segment
requires advancing to its far end. This is a geometric model, not a model of
a particular instrument's cutting mechanism.

Let the line-sphere roots be `s_near <= s_far`, measured from the current tip
along a unit direction. For an outside tip with both roots ahead:

```text
Entry depth = s_near
Core length = s_far - s_near
Required insertion depth = s_far.
```

If the tool has advanced only a distance `d>=0`, the traversed part of that
forward segment is:

```text
max(0, min(d, s_far) - max(0, s_near)), provided s_far >= 0.
```

For an outside tip, this is zero at entry, partial between entry and exit,
and equal to the full available core at `d=s_far`.

As a hand example, choose the tip at `(0,0,0)`, direction `(0,0,1)`, and a
sphere centered at `(0,0,10)` mm with radius 1 mm. Its intersections are at
`s=9` and `s=11`, so the core length is **2 mm** and insertion depth is
**11 mm**. The check `entry + core = depth` gives `9+2=11`.

This is consistent with the professor's clarification about harvesting the
full available core. The model does not justify advancing farther into
surrounding anatomy; it only returns the specified geometric distance.

A tangent has equal roots, so it gives zero core length and a contact depth
if the contact lies ahead. A miss has no insertion depth. A tumor entirely
behind the tip also has no reachable forward core.

## 4. Forward insertion and edge cases

**Student:** How should it handle a tumor behind the tip, a tip inside the tumor, or a reversed insertion direction?

**Chatbot:** Separate the infinite line used for targeting error from the
forward ray used for insertion. Only distances `s>=0` are reachable by
advancing from the current tip.

First normalize the direction and find the ordered signed sphere
intersections. For a real intersection:

| Condition | Meaning | Core length | Insertion depth |
| --- | --- | --- | --- |
| s_far < 0 | Entire sphere intersection is behind the tip | 0 | None |
| s_near < 0 < s_far | Tip is inside; count only tissue ahead | s_far | s_far |
| 0 <= s_near <= s_far | Intersection lies ahead | s_far - s_near | s_far |
| s_near = s_far >= 0 | Forward tangent contact | 0 | Contact distance |

If there is no real intersection, return zero core and `None` for depth.
Do not turn negative depths into positive ones with an absolute value.

The general forward-segment calculation for a reachable intersection is:

```text
entry = max(0, s_near)
exit = max(0, s_far)
core = max(0, exit-entry)
insertion depth = exit.
```

If the tip begins on the near surface and points inward, it can traverse the
forward chord. If it begins on the far surface and points outward, no tissue
remains ahead: both core and depth are zero. This outward-surface case is
not necessarily a tangent; the supporting line can still cross the sphere.

Reversing the direction replaces the roots by `-s_far` and `-s_near`. A sphere
that was entirely ahead becomes entirely behind. The perpendicular distance
to the infinite line stays the same, but the forward plan changes. Do not
silently flip the direction to obtain a successful plan: it may represent
an intentional pose or an upstream axis-sign error that needs attention.

Here are simple hand checks, all with radius 1 mm and tip `(0,0,0)`:

| Sphere center (mm) | Direction | Signed roots (mm) | Expected result |
| --- | --- | --- | --- |
| (0,0,-10) | (0,0,1) | -11, -9 | Behind; core 0, depth None |
| (0,0,0.5) | (0,0,1) | -0.5, 1.5 | Inside; core 1.5 mm, depth 1.5 mm |
| (0,0,10) | (0,0,1) | 9, 11 | Core 2 mm, depth 11 mm |
| (0,0,10) | (0,0,-1) | -11, -9 | Reversed direction; behind |

The implementation uses statuses `intersects`, `tangent`, `miss`, and
`behind`; an inside tip is an intersecting case rather than a separate
`inside` status. The code also uses a small scale-dependent roundoff tolerance
near tangency and at the starting surface. A tolerance should absorb numerical
roundoff without turning a real miss into a hit or removing a small real core.

These calculations are implemented in [planning.py](planning.py). The
required hit, tangent, and miss examples, along with their printed results,
are in the [Q3 README](README.md)
