function [e, c_L, c_Di] = PLLT(b, a0_t, a0_r, c_t, c_r, aero_t, aero_r, geo_t, geo_r, N)
%NewFunction Summary of this function goes here
% Detailed explanation goes here
% Outputs
%   e: span efficiency factor
%   c_L: coefficient of lift
%   c_Di: induced coefficient of drag
% Inputs:
%   b: span [ft]
%   a0_t: cross-sectional lift slope at tip [/rad]
%   a0_r: corss-sectional lift slope at root [/rad]
%   c_t: chord at tip [ft]
%   c_r: chord at root [ft]
%   aero_t: zero lift-angle of attack at tip [deg]
%   aero_r: zero-lift angle of attack at root [deg]
%   geo_t: geometric angle of attack at tip [deg]
%   geo_r: geometric angle of attack at root [deg]
%   N: number of odd terms to include in series expansion
%
% Author: Katherine Korobov
% Collaborators: 
% Date: April 4, 2026

    thetas = linspace(0.0000000000001, pi/2 - 0.000000000001, N);

    geo_aoa_vals = varyGeoAoA(thetas, geo_t, geo_r);
    a0_vals = varyCrossSectionalLiftSlope(thetas, a0_t, a0_r);
    zero_lift_aoa_vals = varyZeroLiftAoA(thetas, aero_t, aero_r);

    LHS = geo_aoa_vals - zero_lift_aoa_vals;

    % build coefficient matrix
    coef = zeros(N, N);

    for i = 1:N
        theta = thetas(i);
        c = varyChordDistribution(theta);

        a0 = a0_vals(i);

        for j = 1:N
            
            m = 2*j - 1; % odd mode index
            
            coef(i,j) = sin(m * theta) * (((2*b)/(a0*c)) + (m / sin(theta)));
            
        end
    end

    A = coef \ LHS;

    % find delta
    sum_term = 0;

    for n = 2:N
        m = 2*n - 1;
        sum_term = sum_term + m * (A(n)/A(1))^2;
    end
    
    AR = ((c_r + c_t) / 2) * b; % calculate aspect ratio assuming linear taper

    e = 1 / (1 + sum_term); % calculate span efficiency factor
    c_L = A(1) * pi * AR;
    c_Di = c_L^2 / (pi * e * AR); 

end


function c = varyChordDistribution(theta, c_r, c_t) 

    if theta >= 0 && theta <= pi/2
        c = c_r - (c_r - c_t) * cos(theta); 
    else 
        c = c_r - (c_t - c_r) * cos(theta);
    end
end

function a0 = varyCrossSectionalLiftSlope(thetas, a0_t, a0_r)
 % * assume linear variation of cross sectional lift slope
    if a0_t == a0_r
        a0 = a0_t * ones(length(thetas), 1);
    else
        a0 = linspace(a0_t, a0_r, length(thetas));
    end
    
end

function aero_aoa = varyZeroLiftAoA(thetas, aero_t, aero_r)
% * assume linear variation of zero lift a0a

    if aero_t == aero_r
        aero_aoa = aero_t * ones(length(thetas), 1);
    else
        aero_aoa = linspace(aero_t, aero_r, length(thetas));
    end
end

function geo_aoa = varyGeoAoA(thetas, geo_t, geo_r)
% ** assume linear variation of geometric angle of attack

    if geo_t == geo_r
        geo_aoa = geo_t * ones(length(thetas),1);

    else 
        geo_aoa = linspace(geo_t, geo_r, length(thetas));
    end
end
