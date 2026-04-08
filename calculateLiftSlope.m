function lift_slope = calculateLiftSlope(alpha, cl)
% calculateLiftSlope finds the cl/alpha lift slope given cl values and
% alpha values
% 
% Assumes  that the lift slope is linear
%
% Author: Katherine Korobov
% Collaborators: 
% Date: 4/8/2026

    p = polyfit(alpha, cl, 1);
    lift_slope = p(1);
    lift_slope = lift_slope * (180/pi); % convert to /rad

end