function thin_airfoil_lift_cl = calculateThinAirfoilCL(alpha)
% calculateThinAirfoilCL finds the sectional lift coefficient from Thin
% Airfoil Theory
% 
%
% Author: Katherine Korobov
% Collaborators: 
% Date: 4/8/2026
    % alpha submited in degrees
    alpha = deg2rad(alpha); % Convert alpha from degrees to radians
    thin_airfoil_lift_cl = (2 * pi) * alpha;
end