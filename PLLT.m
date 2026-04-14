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
% Collaborators: Corey Hannum, John Heflin, Kiana Watson
% Date: April 4, 2026

    thetas = linspace(0.0000000000001, (pi/2) - 0.000000000001, N);

    % build coefficient matrix
    coef = zeros(N, N);

    for i = 1:N
        theta = thetas(i);
        c = varyChordDistribution(theta, c_r, c_t);
        geo_aoa = varyGeoAoA(theta, geo_t, geo_r);
        a0 = varyCrossSectionalLiftSlope(theta, a0_t, a0_r);
        zero_lift_aoa = varyZeroLiftAoA(theta, aero_t, aero_r);

        LHS(i) = geo_aoa - zero_lift_aoa;

        for j = 1:N
            
            m = 2*j - 1; % odd mode index
            
            coef(i,j) = sin(m * theta) * (((2*b)/(a0*c)) + (m / sin(theta)));
            
        end
    end

    A = coef \ LHS';

    % find delta
    sum_term = 0;

    for n = 2:N
        m = 2*n - 1;
        sum_term = sum_term + m * (A(n)/A(1))^2;
    end
    
    S = (b/2) * (c_r + c_t);
    AR = b^2 / S;

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

function a0 = varyCrossSectionalLiftSlope(theta, a0_t, a0_r)
 % * assume linear variation of cross sectional lift slope
    if theta >= 0 && theta <= pi/2
        a0 = a0_r - (a0_r - a0_t) * cos(theta); 
    else 
        a0 = a0_r - (a0_t - a0_r) * cos(theta);
    end
    
end

function aero_aoa = varyZeroLiftAoA(theta, aero_t, aero_r)
% * assume linear variation of zero lift a0a

    if theta >= 0 && theta <= pi/2
        aero_aoa = aero_r - (aero_r - aero_t) * cos(theta); 
    else 
        aero_aoa = aero_r - (aero_t - aero_r) * cos(theta);
    end
end

function geo_aoa = varyGeoAoA(theta, geo_t, geo_r)
% ** assume linear variation of geometric angle of attack
    if theta >= 0 && theta <= pi/2
        geo_aoa = geo_r - (geo_r - geo_t) * cos(theta); 
    else 
        geo_aoa = geo_r - (geo_t - geo_r) * cos(theta);
    end
end
