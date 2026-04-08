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
alpha_zero_lift = 0;
v_inf = 50;
c_l_actual = 2*pi*((12*pi/180)-alpha_zero_lift) % C_l = 2pi*(alpha-(zero lift AoA))
N = 2;
[x_b,y_b] = NACA_airfoils(m,p,t,c,N);

if mod(N,2)==0
    xb = [flip(x_b(N+1:2*N)); x_b(2:N-1,1)];
    yb = [flip(y_b(N+1:2*N)); y_b(2:N-1,1)];
else
    xb = [flip(x_b(N+2:2*N)); x_b(1:N-1,1)];
    yb = [flip(y_b(N+2:2*N)); y_b(1:N-1,1)];
end

CL(N) = Vortex_Panel(xb,yb,v_inf,alpha);
error(N) = (CL(N) - c_l_actual)/c_l_actual;

while abs(error(N)) > 0.01
    N = N + 1;
    [x_b,y_b] = NACA_airfoils(m,p,t,c,N);
    if mod(N,2)==0
        xb = [flip(x_b(N+1:2*N)); x_b(2:N-1,1)];
        yb = [flip(y_b(N+1:2*N)); y_b(2:N-1,1)];
    else
        xb = [flip(x_b(N+2:2*N)); x_b(1:N-1,1)];
        yb = [flip(y_b(N+2:2*N)); y_b(1:N-1,1)];
    end
    CL(N) = Vortex_Panel(xb,yb,v_inf,alpha);
    error(N) = (CL(N) - c_l_actual)/c_l_actual;
end

figure()
plot(2:N,error(2:N)*100, LineWidth=2)
hold on
grid on
yline(-1, LineWidth=1.5)
legend('Error', '1% error')
xlabel('Number of Panels')
ylabel('Percent Error')
ylim([-250, 50])
title('Convergence of the predicted sectional coefficient of lift (c_l) with respect to number of panels (N)')


%% Task 3

param_0012 = struct("m", 0, "p", 0, "t", 12/100);
param_0006 = struct("m", 0, "p", 0, "t", 6/100);
param_0018 = struct("m", 0, "p", 0, "t", 18/100);

num_panels = 50; % will change after task 2 finished

% build airfoils
[x_0006,y_0006] = NACA_airfoils(param_0006.m, param_0006.p, param_0006.t, c, num_panels);
[x_0012,y_0012] = NACA_airfoils(param_0012.m, param_0012.p, param_0012.t, c, num_panels);
[x_0018,y_0018] = NACA_airfoils(param_0018.m, param_0018.p, param_0018.t, c, num_panels);

% we need to flip the x and y to go from trailing edge to leading edge
% clockwise
[flipped_x_0006, fipped_y_0006] = flipPositions(x_0006, y_0006);
[flipped_x_0012, fipped_y_0012] = flipPositions(x_0012, y_0012);
[flipped_x_0018, fipped_y_0018] = flipPositions(x_0018, y_0018);

% collect experimental data
data_0006 = load("NACA_0006_data.mat");
data_0012 = load("NACA_0012_data.mat");

experimental_0006 = data_0006.data;
experimental_0012 = data_0012.data;

% create range for alpha
alpha = linspace(-10, 10, 20); % AoA from 0 deg to 45 deg

for i = 1:length(alpha)
      
    cl_0006(i) = Vortex_Panel(x_0006, y_0006, 1, alpha(i));
    cl_0012(i) = Vortex_Panel(x_0012, y_0012, 1, alpha(i));
    cl_0018(i) = Vortex_Panel(x_0018, y_0018, 1, alpha(i));

end

figure();
hold on;
plot(alpha, cl_0006);
plot(alpha, cl_0012);
plot(alpha, cl_0018);
plot(experimental_0006(:, 1), experimental_0006(:, 2));
plot(experimental_0012(:, 1), experimental_0012(:, 2));
hold off;
legend("Predicted NACA 0006", "Predicted NACA 0012", "Predicted NACA");