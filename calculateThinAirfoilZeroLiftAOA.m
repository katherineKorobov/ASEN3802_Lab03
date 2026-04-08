function zero_lift_alpha = calculateThinAirfoilZeroLiftAOA(airfoil_param, c)
% calculateThinAirfoilZeroLiftAOA determines the angle of attack at L = 0
% using the Thin Airfoil Theory prediction.
% 
% Author: Katherine Korobov
% Collaborators: Kiana Watson
% Date: 4/8/2026

    m = airfoil_param.m;
    p = airfoil_param.p;

    theta = linspace(0, pi, 200);

% x location along chord
x = (c/2) * (1 - cos(theta));
xc = x / c;   % nondimensional x/c

% camber-line slope dz/dx
dzdx = zeros(size(xc));

for i = 1:length(xc)
    if xc(i) < p
        dzdx(i) = (2*m/p^2) * (p - xc(i));
    else
        dzdx(i) = (2*m/(1-p)^2) * (p - xc(i));
    end
end

% Thin airfoil theory zero-lift AoA in radians
zero_lift_alpha = -(1/pi) * trapz(theta, dzdx .* (1 - cos(theta)));

end
