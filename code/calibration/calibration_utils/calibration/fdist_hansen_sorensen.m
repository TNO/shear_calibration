% Hansen, P. F., & Sørensen, J. D. (2002). Reliability-based code 
% calibration of partial safety factors. Paper presented at the Workshop 
% on Reliability Based Code Calibration, Zurich, Switzerland
function dist = fdist_hansen_sorensen(x, x_target)
    dx = x - x_target;
    dist = 4.35 * dx + exp(-4.35 * dx) - 1;
end