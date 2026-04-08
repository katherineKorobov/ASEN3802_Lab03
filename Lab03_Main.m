clc; clear; close all;
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
param_0012 = struct("m", 0, "p", 0, "t", 0.12 * c);

[x_b, y_b] = NACA_airfoils(param_0012.m, param_0012.p, param_0012.t, c, N); 
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
param_2421 = struct("m",0.02 * c, "p", 0.40 * c ,"t", 0.21 * c); 

[x_b2, y_b2] = NACA_airfoils(param_2421.m, param_2421.p, param_2421.t, c, N); 
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
alpha = 12; % degrees
alpha_zero_lift = 0; % 0 for symmetric airfoils
v_inf = 50;
c_l_actual = 2*pi*((12*pi/180)-alpha_zero_lift) % C_l = 2pi*(alpha-(zero lift AoA))
N = 2;
[x_b,y_b] = NACA_airfoils(param_0012.m, param_0012.p, param_0012.t,c,N);

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
    [x_b,y_b] = NACA_airfoils(param_0012.m, param_0012.p, param_0012.t,c,N);
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
plot((2:N).*2,error(2:N)*100, LineWidth=2)
hold on
grid on
yline(-1, LineWidth=1.5)
legend('Error', '1% error')
xlabel('Number of Panels')
ylabel('Percent Error')
ylim([-250, 50])
title('Convergence of the predicted sectional coefficient of lift (c_l) with respect to number of panels (N)')

%% Task 3
% create range for alpha
alpha = linspace(-10, 10, 20); % [deg]

% Vortex Panel Method
% * NACA 0012 is defined above earlier
param_0006 = struct("m", 0, "p", 0, "t", 0.06 * c);
param_0018 = struct("m", 0, "p", 0, "t", 0.18 * c);

num_panels = 28; % from Task 2

% build airfoils
[x_0006, y_0006] = NACA_airfoils(param_0006.m, param_0006.p, param_0006.t, c, num_panels);
[x_0012, y_0012] = NACA_airfoils(param_0012.m, param_0012.p, param_0012.t, c, num_panels);
[x_0018, y_0018] = NACA_airfoils(param_0018.m, param_0018.p, param_0018.t, c, num_panels);

% we need to flip the x and y to go from trailing edge to leading edge
% clockwise
[flipped_x_0006, flipped_y_0006] = flipPositions(x_0006, y_0006);
[flipped_x_0012, flipped_y_0012] = flipPositions(x_0012, y_0012);
[flipped_x_0018, flipped_y_0018] = flipPositions(x_0018, y_0018);

for i = 1:length(alpha)     
    cl_0006(i) = Vortex_Panel(flipped_x_0006, flipped_y_0006, 1, alpha(i));
    cl_0012(i) = Vortex_Panel(flipped_x_0012,flipped_y_0012, 1, alpha(i));
    cl_0018(i) = Vortex_Panel(flipped_x_0018, flipped_y_0018, 1, alpha(i));
end

zero_lift_aoa_0006 = interp1(cl_0006, alpha, 0, "linear"); % find alpha for cl = 0 -> L = 0
zero_lift_aoa_0012 = interp1(cl_0012, alpha, 0, "linear");
zero_lift_0018 = interp1(cl_0018, alpha, 0, "linear");

lift_slope_0006 = calculateLiftSlope(alpha, cl_0006);
lift_slope_0012 = calculateLiftSlope(alpha, cl_0012);
lift_slope_0018 = calculateLiftSlope(alpha, cl_0018);

% Experimental Method
data_0006 = load("NACA_0006_data.mat"); % collect experimental data
data_0012 = load("NACA_0012_data.mat"); % collect experimental data
experimental_0006 = data_0006.data;
experimental_0012 = data_0012.data;

zero_lift_aoa_exp_0006 = interp1(experimental_0006(:, 2), experimental_0006(:,1), 0, "linear");
zero_lift_aoa_exp_0012 = interp1( experimental_0012(:, 2), experimental_0012(:,1), 0, "linear");

lift_slope_exp_0006 = calculateLiftSlope(experimental_0006(:, 1), experimental_0006(:, 2));
lift_slope_exp_0012 = calculateLiftSlope(experimental_0012(:, 1), experimental_0012(:, 2));

% Thin Airfoil Theory (TAT)
TAT_cl_0006 = calculateThinAirfoilCL(alpha, param_0006, c);
TAT_cl_0012 = calculateThinAirfoilCL(alpha, param_0012, c);
TAT_cl_0018 = calculateThinAirfoilCL(alpha, param_0018, c);

TAT_zero_lift_aoa_0006 = calculateThinAirfoilZeroLiftAOA(param_0006, c);
TAT_zero_lift_aoa_0012 = calculateThinAirfoilZeroLiftAOA(param_0012, c);
TAT_zero_lift_aoa_0018 = calculateThinAirfoilZeroLiftAOA(param_0018, c);

[TAT_lift_slope_0006, TAT_lift_slope_0012, TAT_lift_slope_0018] = deal((2*pi) * (pi / 180)) ;

% Combine
figure();
hold on;
plot(alpha, cl_0006,"-" ,"LineWidth", 1.5);
plot(alpha, cl_0012, "-" ,"LineWidth", 1.5);
plot(alpha, cl_0018, "-","LineWidth", 1.5);
plot(experimental_0006(:, 1), experimental_0006(:, 2), ":" ,"LineWidth", 1.5);
plot(experimental_0012(:, 1), experimental_0012(:, 2), ":", "LineWidth", 1.5);
plot(alpha, TAT_cl_0006, "--", "LineWidth", 1.5);
plot(alpha, TAT_cl_0012, "--", "LineWidth", 1.5);
plot(alpha, TAT_cl_0018, "--", "LineWidth", 1.5);
hold off;
legend("Vortex Panel NACA 0006", "Vortex Panel NACA 0012", "Vortex Panel NACA 0018", ...
    "Experimental NACA 0006", "Experimental NACA 0012", "TAT NACA 0006", "TAT NACA 0012", "TAT NACA 0018", "Location","southeast");
xlim([-11, 11]);
xlabel("Angle of Attack [deg]");
ylabel("Sectional Coefficient of Lift");
title("Sectional Coefficient of Lift v. Angle of Attack for Varying Airfoil Data");


%% Task 4: Effect of Airfoil Camber on Lift

% create range for alpha
alpha = linspace(-10, 10, 20); % [deg]

num_panels = 28; % from Task 2

% Airfoil params
param_0012 = struct("m", 0, "p", 0, "t", 0.12 * c);
param_2412 = struct("m", 0.02, "p", 0.4, "t", 0.12 * c);
param_4412 = struct("m", 0.04, "p", 0.4, "t", 0.12 * c);

% Build airfoils
[x_0012, y_0012] = NACA_airfoils(param_0012.m, param_0012.p, param_0012.t, c, num_panels);
[x_2412, y_2412] = NACA_airfoils(param_2412.m, param_2412.p, param_2412.t, c, num_panels);
[x_4412, y_4412] = NACA_airfoils(param_4412.m, param_4412.p, param_4412.t, c, num_panels);

% Reorder points clockwise 
[flipped_x_0012, flipped_y_0012] = flipPositions(x_0012, y_0012);
[flipped_x_2412, flipped_y_2412] = flipPositions(x_2412, y_2412);
[flipped_x_4412, flipped_y_4412] = flipPositions(x_4412, y_4412);

% Vortex Panel
cl_0012 = zeros(1, length(alpha));
cl_2412 = zeros(1, length(alpha));
cl_4412 = zeros(1, length(alpha));

for i = 1:length(alpha)
    cl_0012(i) = Vortex_Panel(flipped_x_0012, flipped_y_0012, 1, alpha(i));
    cl_2412(i) = Vortex_Panel(flipped_x_2412, flipped_y_2412, 1, alpha(i));
    cl_4412(i) = Vortex_Panel(flipped_x_4412, flipped_y_4412, 1, alpha(i));
end

% Zero-Lift AoA for VP 
zero_lift_aoa_0012 = interp1(cl_0012, alpha, 0, "linear");
zero_lift_aoa_2412 = interp1(cl_2412, alpha, 0, "linear");
zero_lift_aoa_4412 = interp1(cl_4412, alpha, 0, "linear");

% Lift slope for VP
lift_slope_0012 = calculateLiftSlope(alpha, cl_0012);
lift_slope_2412 = calculateLiftSlope(alpha, cl_2412);
lift_slope_4412 = calculateLiftSlope(alpha, cl_4412);

% Experimental Data
data_0012 = load("NACA_0012_data.mat");
data_2412 = load("NACA_2412_data.mat");
data_4412 = load("NACA_4412_data.mat");

experimental_0012 = data_0012.data;
experimental_2412 = data_2412.data;
experimental_4412 = data_4412.data;

% Zero-lift AoA exp
zero_lift_aoa_exp_0012 = interp1(experimental_0012(:,2), experimental_0012(:,1), 0, "linear");
zero_lift_aoa_exp_2412 = interp1(experimental_2412(:,2), experimental_2412(:,1), 0, "linear");
zero_lift_aoa_exp_4412 = interp1(experimental_4412(:,2), experimental_4412(:,1), 0, "linear");

% Lift slope exp
lift_slope_exp_0012 = calculateLiftSlope(experimental_0012(:,1), experimental_0012(:,2));
lift_slope_exp_2412 = calculateLiftSlope(experimental_2412(:,1), experimental_2412(:,2));
lift_slope_exp_4412 = calculateLiftSlope(experimental_4412(:,1), experimental_4412(:,2));

%Thin airfoil theory (TAT)
%alpha_L0_0012 = calculateThinAirfoilZeroLiftAOA(param_0012, c);
zero_lift_aoa_TAT_0012 = 0; 
zero_lift_aoa_TAT_2412 = calculateThinAirfoilZeroLiftAOA(param_2412, c);
zero_lift_aoa_TAT_4412 = calculateThinAirfoilZeroLiftAOA(param_4412, c);   

alpha_rad = deg2rad(alpha);

cl_TAT_0012 = 2*pi*(alpha_rad - zero_lift_aoa_TAT_0012);
cl_TAT_2412 = 2*pi*(alpha_rad - zero_lift_aoa_TAT_2412);
cl_TAT_4412 = 2*pi*(alpha_rad - zero_lift_aoa_TAT_4412);

lift_slope_TAT = 2*pi/180; 

% Combined Plot
figure();
hold on;
plot(alpha, cl_0012, 'b:', "LineWidth", 2);
plot(alpha, cl_2412, 'g:', "LineWidth", 2);
plot(alpha, cl_4412, 'm:', "LineWidth", 2);
plot(alpha, cl_TAT_0012, "b", "LineWidth", 1.5);
plot(alpha, cl_TAT_2412, "g", "LineWidth", 1.5);
plot(alpha, cl_TAT_4412, "m", "LineWidth", 1.5);
plot(experimental_0012(:,1), experimental_0012(:,2), 'b--', "LineWidth", 1.5);
plot(experimental_2412(:,1), experimental_2412(:,2), 'g--', "LineWidth", 1.5);
plot(experimental_4412(:,1), experimental_4412(:,2), 'm--', "LineWidth", 1.5);
hold off;
legend("Predicted NACA 0012", "Predicted NACA 2412", "Predicted NACA 4412", ...
    "TAT NACA 0012", "TAT NACA 2412", "TAT NACA 4412", ...
    "Experimental NACA 0012", "Experimental NACA 2412", "Experimental NACA 4412", ...
    "Location", "southeast");
xlim([-11, 11]);
xlabel("Angle of Attack [deg]");
ylabel("Sectional Coefficient of Lift");
title("Sectional Coefficient of Lift v. Angle of Attack for Varying Airfoil Camber");


% Zero-Lift AoA table
fprintf('\n--- Zero-Lift Angle of Attack [deg] ---\n');
fprintf('%-12s %-12s %-12s %-12s\n', 'Airfoil', 'VPM', 'TAT', 'Experimental');
fprintf('%-12s %-12.4f %-12.4f %-12.4f\n', 'NACA 0012', zero_lift_aoa_0012, zero_lift_aoa_TAT_0012, zero_lift_aoa_exp_0012);
fprintf('%-12s %-12.4f %-12.4f %-12.4f\n', 'NACA 2412', zero_lift_aoa_2412, zero_lift_aoa_TAT_2412, zero_lift_aoa_exp_2412);
fprintf('%-12s %-12.4f %-12.4f %-12.4f\n', 'NACA 4412', zero_lift_aoa_4412, zero_lift_aoa_TAT_4412, zero_lift_aoa_exp_4412);

% Lift Slope table
fprintf('\n--- Lift Slope [per deg] ---\n');
fprintf('%-12s %-12s %-12s %-12s\n', 'Airfoil', 'VPM', 'TAT', 'Experimental');
fprintf('%-12s %-12.4f %-12.4f %-12.4f\n', 'NACA 0012', lift_slope_0012, lift_slope_TAT, lift_slope_exp_0012);
fprintf('%-12s %-12.4f %-12.4f %-12.4f\n', 'NACA 2412', lift_slope_2412, lift_slope_TAT, lift_slope_exp_2412);
fprintf('%-12s %-12.4f %-12.4f %-12.4f\n', 'NACA 4412', lift_slope_4412, lift_slope_TAT, lift_slope_exp_4412);
