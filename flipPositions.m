function [x_flipped,y_flipped] = flipPositions(x,y)
% flipPositions flips the positions of the x and y coordinates for Vortex Panel Method.
% 
% Vortext panel method requires the x and y location sto be arranged in a 
% clockwise order starting from teh trailing edge
%
% Author: Katherine Korobov
% Collaborators: 
% Date: 4/8/2026


N = length(x)/2;

if mod(N,2) == 0
    x_flipped = [flip(x(N+1:2*N)); x(2:N-1,1)];
    y_flipped = [flip(y(N+1:2*N)); y(2:N-1,1)];
else
    x_flipped = [flip(x(N+2:2*N)); x(1:N-1,1)];
    y_flipped = [flip(y(N+2:2*N)); y(1:N-1,1)];
end



end