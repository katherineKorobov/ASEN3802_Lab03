function [e, c_L, c_Di] = PLLT(b, a0_t, a0_r, c_t, c_r, aero_t, aero_r, geo_t, geo_r, N)
% Provides span efficiency, coefficient of lift, and induced coefficient of
% drag using Prandtl Lifting Line Theory for varying types of 3D wings.
% 
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
% Date: April 14, 2026
    
    aero_t = deg2rad(aero_t); % convert to radians
    aero_r = deg2rad(aero_r); % convert to radians
    geo_t  = deg2rad(geo_t); % convert to radians
    geo_r  = deg2rad(geo_r); % convert to radians

    S = (b/2) * (c_r + c_t); % compute surface area
    AR = b^2 / S; % compute aspect ratio
    
    thetas = (1:N) * pi / (2*N); % build range of theta's across wing
    coef = zeros(N, N); % build coefficient matrix

    for i = 1:N
        theta = thetas(i);
        c = varyChordDistribution(theta, c_r, c_t);
        geo_aoa = varyGeoAoA(theta, geo_t, geo_r);
        a0 = varyCrossSectionalLiftSlope(theta, a0_t, a0_r);
        zero_lift_aoa = varyZeroLiftAoA(theta, aero_t, aero_r);

        LHS(i) = geo_aoa - zero_lift_aoa;

        for j = 1:N
            m = 2*j - 1; % odd mode index
            coef(i,j) = sin(m * theta) * (((4*b)/(a0*c)) + (m / sin(theta)));     
        end
    end

    A = coef \ LHS';

    
    sum_term = 0;
    for n = 2:N
        m = 2*n - 1;
        sum_term = sum_term + m * (A(n)/A(1))^2; % find delta
    end

    e = 1 / (1 + sum_term); % calculate span efficiency factor
    c_L = A(1) * pi * AR; % compute sectional lift
    c_Di = c_L^2 / (pi * e * AR); % compute iduced drag coefficient

end

function c = varyChordDistribution(theta, c_r, c_t) 
% Calculates the chord length for a given theta
% Outputs:
%   c : chord length [ft]
% Inputs:
%   theta: position on wing [rad]
%   c_r: root chord length [ft]
%   c_t: tip chord length [ft]
%
% Author: Katherine Korobov
% Collaborators: Corey Hannum, John Heflin, Kiana Watson
% Date: April 14, 2026

    c = c_r - (c_r - c_t) * abs(cos(theta)); 
end

function a0 = varyCrossSectionalLiftSlope(theta, a0_t, a0_r)
% Calculates the chord length for a given theta
% Outputs:
%   a0 : cross-sectional lift slope [/rad]
% Inputs:
%   theta: position on wing [rad]
%   a0_t: root cross-sectional lift slope [/rad]
%   a0_r: tip cross-sectional lift slope [/rad]
%
% Assumes linear variation across wing with chord
%
% Author: Katherine Korobov
% Collaborators: Corey Hannum, John Heflin, Kiana Watson
% Date: April 14, 2026

    a0 = a0_r - (a0_r - a0_t) * abs(cos(theta)); 
end

function aero_aoa = varyZeroLiftAoA(theta, aero_t, aero_r)
% Calculates the zero-lift angle of attack for a given theta
% Outputs:
%   aero_aoa : zero-lift angle of attack [rad]
% Inputs:
%   theta: position on wing [rad]
%   aero_t: root zero-lift angle of attack [rad]
%   aero_r: tip zero-lift angle of attack [rad]
%
% Assumes linear variation across wing with chord
%
% Author: Katherine Korobov
% Collaborators: Corey Hannum, John Heflin, Kiana Watson
% Date: April 14, 2026

    aero_aoa = aero_r - (aero_r - aero_t) * abs(cos(theta)); 
end

function geo_aoa = varyGeoAoA(theta, geo_t, geo_r)
% Calculates the geometric angle of attack for a given theta
% Outputs:
%   geo_aoa : zero lift angle of attack [rad]
% Inputs:
%   theta: position on wing [rad]
%   geo_t: root geometric angle of attack[rad]
%   geo_r: tip geometric angle of attack [rad]
%
% Assumes linear variation across wing with chord
%
% Author: Katherine Korobov
% Collaborators: Corey Hannum, John Heflin, Kiana Watson
% Date: April 14, 2026

    geo_aoa = geo_r - (geo_r - geo_t) * abs(cos(theta)); 
end
