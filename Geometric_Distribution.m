function x=Geometric_Distribution(c,n)
% Inputs: chord length and number of points (c,n)
% Outputs: x (equiangularly distributed
% Function takes the chord length (c) and the number of points needed (n)
% and creates a x distribution using the equiangular rays method about the
% center of the chord (c/2)
%
% Author: {John Heflin}
% Date: {03/31/2026}

theta=linspace(0,pi,n); %radians


r=c./2; %radius=chord/2 (m)

x=r-r.*cos(theta) %x coordinates using equiangular rays (start r distance away and subtract the projection of the ray onto the x axis)(m)
end

%checking distribution
%figure
%plot(x,1,'*')

