function zero_lift_alpha = calculateThinAirfoilZeroLiftAOA(airfoil_param, c)
% calculateThinAirfoilZeroLiftAOA determines the angle of attack at L = 0
% using the Thin Airfoil Theory prediction.
% 
% Author: Katherine Korobov
% Collaborators: 
% Date: 4/8/2026

    theta = linspace(0, pi, 200);

    x = (c / 2) * (1- cos(theta));

    dzdx = zeros(length(x));

    % Calculate the slope of the lift curve
    for i = 1:length(x)
        dzdx(i) =  
    end

    integrand = dzdx .* (cos(theta) - 1);
    zero_lift_alpha = - (1/pi) * trapz(theta, integrand);

end