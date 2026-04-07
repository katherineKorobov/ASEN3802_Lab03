function [x_flipped,y_flipped] = flipPositions(x,y)

% Vortext panel method requires the x and y location sto be arranged in a 
% clockwise order starting from teh trailing edge

trailing_edge_idx = length(x)/2;

x_half = x(1:trailing_edge_idx); % go from beginning ot trailing edge 0 -> 1
y_half = y(1:trailing_edge_idx); % corresponding y values up to the trailing edge

x_flipped = flip(x_half); % flip the x positions
y_flipped = flip(y_half); % flip the corresponding y positions

x_flipped = [-x_flipped, x_half];
y_flipped = [-y_flipped, y_half];



end