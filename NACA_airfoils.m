function [x_b, y_b] = NACA_airfoils(m,p,t,c,N)
    
%NACA_airfoils outputs the coordinates of boundary points for a NACA airfoil.  
    % The function takes in the max camber, the location of max camber, the
    % thickness, the chord length, and the number of employed panels. It
    % then calculates the x and y coordinates of the upper and lower
    % airfoil surfaces. It then outputs these in the vectors x_b and y_b,
    % which represent the x-locations and y-locations of the boundary
    % points respectively. 
    %
% Author: Kiana Watson
% Collaborators: 
% Date: {should include the date last revised}

x_u = zeros(N, 1); % n is total number of x's being used
x_l = zeros(N, 1);
y_u = zeros(N, 1); 
y_l = zeros(N, 1);

x=Geometric_Distribution(c,N);
for i = 1 : N %adjust how x is obtained based on method used 
    

    x=Geometric_Distribution(c,N);
    
    y_t = (t/0.2).*c .* ((0.2969.*sqrt(x(i)./c)) - (0.1260.*(x(i)./c)) - (0.3516.*(x(i)./c).^2) + (0.2843.*(x(i)./c).^3) - (0.1036.*(x(i)./c).^4));

    if (0 <= x(i)) && (x(i) < (p.*c))

        y_c = m.*(x(i)/(p.^2)) .* (2.*p - (x(i)./c)); 
        dy_c = (m.*(2.*p - x(i)./c))/p.^2 - (m.*x(i))./(c.*p.^2); %derivative with respect to x 
        
    elseif ((p.*c)<= x(i)) && (x(i) <= c)

        y_c = m.*((c-x(i))./(1-p).^2) .* (1 + (x(i)./c) - (2.*p)); 
        dy_c = (m.*(c - x(i)))./(c.*(p - 1).^2) - (m.*(x(i)./c - 2.*p + 1))/(p - 1).^2; 

    else
        fprintf("x not within bounds") 
    end 
    
    xi = atan(dy_c); 

    x_u(i) = x(i) - y_t.*sin(xi);
    x_l(i) = x(i) + y_t.*sin(xi); 
    y_u(i) = y_c + y_t.*cos(xi); 
    y_l(i) = y_c - y_t.*cos(xi);

end 

x_b = [x_u; x_l]; 
y_b = [y_u; y_l]; 
    
end

% Note: coordinates must go clockwise from tailing edge 
