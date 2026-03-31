function x=Geometric_Distribution(c,n)
theta=linspace(0,pi,n); %radians


r=c./2; %radius=chord/2 (m)

x=r-r.*cos(theta) %x coordinates using equiangular rays (start r distance away and subtract the projection of the ray onto the x axis)(m)
end

%checking distribution
%figure
%plot(x,1,'*')

