# CT-guided Bone Biopsy Planning

CISC/CMPE 330 Assignment 2: calibration, registration, geometric planning, and targeting-error analysis for a tracked biopsy needle.

## Project components

| Component | Purpose | Repository status |
| --- | --- | --- |
| Tool calibration | Estimate the tool-frame tip and axis from constrained tracked movements; evaluate exact and noisy observations. | Separate team contribution; not included in this version. |
| Registration | Recover the tool-to-CT rigid transformation and transform the tip and axis. | Implemented and tested. |
| Biopsy planning | Compute targeting error, available tissue-core length, and insertion depth for a spherical tumor. | Implemented and tested. |
| Targeting-error simulation | Propagate marker-localization and target-displacement errors and estimate tolerable target uncertainty. | Separate team contribution; not included in this version. |

The shared functions support the complete workflow:

```text
Calibration -> tool-frame tip and axis
                            |
Tool markers + CT markers -> Registration -> CT tip and axis -> Biopsy plan
                                   |
                 Noisy markers -> repeated registration -> targeting-error analysis
```


## Coordinate and geometry conventions

- Positions, distances, and radii are in **millimetres**.
- Marker arrays have shape `(N, 3)`, with one point per row. Matching rows must
  identify the same marker in both frames. Registration does not infer matches.
- Tip and axis inputs have shape `(3,)`. Axes are normalized and point in the
  direction of insertion.
- A transform maps **tool/source to CT/destination**:
  `point_ct = R @ point_tool + t` and `axis_ct = R @ axis_tool`.
- Targeting error is the shortest distance from the tumor center to the infinite
  needle axis. Core length and insertion depth use the forward ray from the tip.
- Full-core insertion depth reaches the far tumor boundary. Tangency gives zero
  core and a finite contact depth. A miss or a target wholly behind the tip gives
  zero core and `None` for insertion depth.
- A tip already inside the sphere counts only the tissue remaining ahead.

Figure 2 in the assignment does not label the coordinate units.
`assignment_tool_model` therefore requires `mm_per_figure_unit`. The synthetic
examples choose `1.0`; the physical conversion must be agreed before interpreting
the team's uncertainty results. Direction vectors are not scaled.

## Example

```python
from bone_biopsy import assignment_tool_model, register_tool, plan_biopsy

# Choose an explicit scale for the synthetic example.
model = assignment_tool_model(mm_per_figure_unit=1.0)
ct_markers = model.markers + [10.0, -5.0, 2.0]

# Register the tool using matching marker observations.
tool = register_tool(model.markers, ct_markers, model.tip, model.axis)

# Plan a trajectory through a sphere 100 mm ahead of the tip.
plan = plan_biopsy(tool.tip, tool.axis, [10.0, -105.0, 2.0], 5.0)
print(plan.targeting_error, plan.tissue_core_length, plan.insertion_depth)
# Approximately: 0.0, 10.0, 105.0
```

For team integration, replace the ground-truth tip and axis with calibrated
values expressed in the same model frame as the tool markers.


## Public functions



## View the solution



## Run the implementation




## References