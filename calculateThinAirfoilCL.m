function thin_airfoil_lift_cl = calculateThinAirfoilCL(alpha, airfoil_param, c)
% calculateThinAirfoilCL finds the sectional lift coefficient from Thin
% Airfoil Theory
% 
%
% Author: Katherine Korobov
% Collaborators: 
% Date: 4/8/2026

    % alpha submited in degrees
    alpha_rad = deg2rad(alpha); % Convert alpha from degrees to radians
    alpha_zero_lift = calculateThinAirfoilZeroLiftAOA(airfoil_param, c);
    thin_airfoil_lift_cl = (2 * pi) .* (alpha_rad - alpha_zero_lift);
end