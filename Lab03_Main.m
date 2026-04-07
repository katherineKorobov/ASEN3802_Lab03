clc
clear
close all
%% ASEN 3802 - Lab 03 - Main
% Part 1: Analysis of 2D Airfoils
%   Task 1: NACA 4-digit Airfoil Generator
%   Task 2: Convergence study on number of panels using vortex panel method
%
% Authors: {Corey Hannum, John Heflin, Katherine Korobov, and Kiana Watson}
% Date: {03/31/2026}


%% Part 1




%%Task 1

c=1; % chord length (m)
N=50; % number of panels

%% Task 1
c=1; % chord length (m)
N=50; % number of panels

%NACA 0021 ----------------
m = 0;
p = 0;
t = 0.21 .* c; 

[x_b, y_b] = NACA_airfoils(m, p, t, c, N); 
figure
hold on; 
plot(x_b(1:N), y_b(1:N),'b') % upper surface
plot(x_b(N+1:2.*N), y_b(N+1:2.*N),'b') % lower surface
xlabel("x")
ylabel("y")
xlim([0,c]);
ylim([-c./2,c./2]);
title("NACA 0021")
hold off;

%NACA 2421 ----------------
m_2 = 0.02 * c;
p_2 = 0.40 * c; 
t_2 = 0.21 * c; 

[x_b2, y_b2] = NACA_airfoils(m_2, p_2, t_2, c, N); 
y_c=y_b2(1:N)+y_b2(N+1:2.*N);
figure
hold on; 
plot(x_b2(1:N), y_b2(1:N),'b') % upper surface
plot(x_b2(N+1:2.*N), y_b2(N+1:2.*N),'b') % lower surface
plot(x_b2(1:N),y_c,'k')
xlabel("x")
ylabel("y")
xlim([0,c]);
ylim([-c./2,c./2]);
legend('NACA 2421','NACA 2421','Mean chamberline')
title("NACA 2421")
hold off; 


%% task 2
% NACA 0012
m = 0; 
p = 0;
t = 12/100;
alpha = 12; % degrees
v_inf = 50;
error = 100;
c_l_actual = 1;
N = 1;
[x_b,y_b] = NACA_airfoils(m,p,t,c,N);
CL(N) = Vortex_Panel(x_b,y_b,v_inf,alpha);
error(N) = (CL(N) - c_l_actual)/c_l_actual;

while error(N) > 0.01
    N = N + 1;
    [x_b,y_b] = NACA_airfoils(m,p,t,c,N);
    CL(N) = Vortex_Panel(x_b,y_b,v_inf,alpha);
    error(N) = (CL(N) - c_l_actual)/c_l_actual;
end

figure()
plot(1:N,error)
hold on
yline(0.01)
xlabel('Number of Panels')
ylabel('Percent Error')
title('Convergence of the predicted sectional coefficient of lift (c_l) with respect to number of panels (N)')


