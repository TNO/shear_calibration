# data

* `beam_p_adj.mat` database shared by Yuguang in 2019 October, original name is kept for traceability (intentionally not kept under version control (not shared), because we do not have permission to share it)
* `load_comb_prevalence_weights.xlsx` prevalence weights for all possible load combinations, read and used in calibration

## Shear database

Columns:
* `bw`: [mm] smallest width of the cross section in the tensile area (the web width of the cross section of a beam)
* `d`: [mm] effective depth of the cross section of a beam
* `dg`: [mm] the size parameter describing the failure zone roughness
* `fc`: [MPa] concrete compressive strength
* `mvd`: shear slenderness ratio
* `rho`: [%] longitudinal reinforcement ratio
* `Vu`: [kN] ultimate shear force at the critical cross section where the critical shear crack crosses it neutral axis

## Prevalence weights

* The weights (`w`) on the `traffic` sheet are the same as those reported in the paper (up to a constant scaler that does not affect the optimum location). In the Excel sheet the the values are scaled to reach 1.0 with the largest weight.
* For all load combinations with two loads, e.g. `snow_wind`, the weight matrix is obtained as: `w * w'`.
* The weights are given over an evenly spaced grid of `chi` values: `0:0.1:1`. During the calibration the weights corresponding the `chi` values are computed by linearly interpolating between values over the grid (in 1d and 2d as well).

