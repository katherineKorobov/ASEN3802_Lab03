clc; clear; close all;
%% ASEN 3802 - Lab 03 - Main
% Part 1: Analysis of 2D Airfoils
%   Task 1: NACA 4-digit Airfoil Generator
%   Task 2: Convergence study on number of panels using vortex panel method
%   Task 3: Analysis of the accuracy of thin airfoil theory and vortex
%   panel method for symmetric airfoils
%   Task 4: Analysis of the accuracy of thin airfoil theory and vortex
%   panel method for cambered airfoils
% Authors: {Corey Hannum, John Heflin, Katherine Korobov, and Kiana Watson}
% Date: {04/14/2026}

set(groot, "defaultAxesTickLabelInterpreter","latex"); 
set(groot, "defaultLegendInterpreter","latex");
set(groot, "defaultTextInterpreter", "latex");
set(groot,"defaultAxesfontsize", 14);

%% Part 01
plotNACA0021 = 0; % plot toggle
plotNACA2421 = 0; % plot toggle
plotPredictedCL = 0; % plot toggle
plotThickAirfoil = 0; % plot toggle
plotAirfoilCamber = 0; % plot toggle

%% Task 1
c=1; % chord length (m)
N=50; % number of panels

%NACA 0021 ----------------
param_0012 = struct("m", 0, "p", 0, "t", 0.12 * c);
[x_0021, y_0021] = NACA_airfoils(param_0012.m, param_0012.p, param_0012.t, c, N); 

if plotNACA0021
    figure
    hold on; 
    plot(x_0021(1:N), y_0021(1:N),'b') % upper surface
    plot(x_0021(N+1:2.*N), y_0021(N+1:2.*N),'b') % lower surface
    xlabel("x")
    ylabel("y")
    xlim([0,c]);
    ylim([-c./2,c./2]);
    title("NACA 0021")
    hold off;
end

%NACA 2421 ----------------
param_2421 = struct("m",0.02 * c, "p", 0.40 * c ,"t", 0.21 * c); 
[x_2421, y_2421] = NACA_airfoils(param_2421.m, param_2421.p, param_2421.t, c, N); 
y_c=y_2421(1:N)+y_2421(N+1:2.*N);

if plotNACA2421
    figure
    hold on; 
    plot(x_2421(1:N), y_2421(1:N),'b') % upper surface
    plot(x_2421(N+1:2.*N), y_2421(N+1:2.*N),'b') % lower surface
    plot(x_2421(1:N),y_c,'k')
    xlabel("x")
    ylabel("y")
    xlim([0,c]);
    ylim([-c./2,c./2]);
    legend('NACA 2421','NACA 2421','Mean chamberline')
    title("NACA 2421")
    hold off; 
end

%% task 2
% NACA 0012
alpha = 12; % degrees
alpha_zero_lift = 0; % 0 for symmetric airfoils
v_inf = 50;
c_l_actual = 2*pi*((12*pi/180)-alpha_zero_lift); % C_l = 2pi*(alpha-(zero lift AoA))
N = 2;
[x_0021_2,y_0021_2] = NACA_airfoils(param_0012.m, param_0012.p, param_0012.t,c,N);

% flip the arrays to start at TE and rotate CW
if mod(N,2)==0
    xb = [flip(x_0021_2(N+1:2*N)); x_0021_2(2:N-1,1)];
    yb = [flip(y_0021_2(N+1:2*N)); y_0021_2(2:N-1,1)];
else
    xb = [flip(x_0021_2(N+2:2*N)); x_0021_2(1:N-1,1)];
    yb = [flip(y_0021_2(N+2:2*N)); y_0021_2(1:N-1,1)];
end

CL(N) = Vortex_Panel(xb,yb,v_inf,alpha);
error(N) = (CL(N) - c_l_actual)/c_l_actual;

% loop until error is less than 1%
while abs(error(N)) > 0.01
    N = N + 1;
    [x_0021_2,y_0021_2] = NACA_airfoils(param_0012.m, param_0012.p, param_0012.t,c,N);
    % flip the arrays to start at TE and rotate CW
    if mod(N,2)==0
        xb = [flip(x_0021_2(N+1:2*N)); x_0021_2(2:N-1,1)];
        yb = [flip(y_0021_2(N+1:2*N)); y_0021_2(2:N-1,1)];
    else
        xb = [flip(x_0021_2(N+2:2*N)); x_0021_2(1:N-1,1)];
        yb = [flip(y_0021_2(N+2:2*N)); y_0021_2(1:N-1,1)];
    end
    % calculate Cl and error
    CL(N) = Vortex_Panel(xb,yb,v_inf,alpha);
    error(N) = (CL(N) - c_l_actual)/c_l_actual;
end

if plotPredictedCL
    figure()
    plot((2:N).*2,error(2:N)*100, LineWidth=2)
    hold on
    grid on
    yline(-1, LineWidth=1.5)
    legend('Error', '1% error')
    xlabel('Number of Panels')
    ylabel('Percent Error')
    ylim([-250, 50])
    title('Convergence of the predicted sectional coefficient of lift ($c_l$) with respect to number of panels (N)')
end

%% Task 3
% create range for alpha
alpha = linspace(-10, 10, 200); % [deg]

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
zero_lift_aoa_0012 = interp1(cl_0012, alpha, 0, "linear"); % [deg]
zero_lift_aoa_0018 = interp1(cl_0018, alpha, 0, "linear"); % [deg]

lift_slope_0006 = calculateLiftSlope(alpha, cl_0006);
lift_slope_0012 = calculateLiftSlope(alpha, cl_0012);
lift_slope_0018 = calculateLiftSlope(alpha, cl_0018);

% Experimental Method
data_0006 = load("NACA_0006_data.mat"); % collect experimental data alpha in deg
data_0012 = load("NACA_0012_data.mat"); % collect experimental data alpha in deg
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

TAT_zero_lift_aoa_0006 = rad2deg(calculateThinAirfoilZeroLiftAOA(param_0006, c)); 
TAT_zero_lift_aoa_0012 = rad2deg(calculateThinAirfoilZeroLiftAOA(param_0012, c));
TAT_zero_lift_aoa_0018 = rad2deg(calculateThinAirfoilZeroLiftAOA(param_0018, c));

[TAT_lift_slope_0006, TAT_lift_slope_0012, TAT_lift_slope_0018] = deal(2*pi * (pi/180));

% Combine
if plotThickAirfoil
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
end

%% Task 4: Effect of Airfoil Camber on Lift
% * Uses same alpha and num_panels
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
cl_ave=(cl_0012+cl_2412)./2;
% Zero-Lift AoA for VP 
zero_lift_aoa_0012 = interp1(cl_0012, alpha, 0, "linear"); % [deg]
zero_lift_aoa_2412 = interp1(cl_2412, alpha, 0, "linear"); % [deg]
zero_lift_aoa_4412 = interp1(cl_4412, alpha, 0, "linear"); % [deg]

% Lift slope for VP, a0
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
zero_lift_aoa_TAT_2412 = rad2deg(calculateThinAirfoilZeroLiftAOA(param_2412, c));
zero_lift_aoa_TAT_4412 = rad2deg(calculateThinAirfoilZeroLiftAOA(param_4412, c));   

cl_TAT_0012 = calculateThinAirfoilCL(alpha, param_0012, c);
cl_TAT_2412 = calculateThinAirfoilCL(alpha, param_2412, c);
cl_TAT_4412 = calculateThinAirfoilCL(alpha, param_4412, c);

lift_slope_TAT = 2*pi * (pi / 180); 

% Combined Plot
if plotAirfoilCamber
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
end


%% Part 02
plotInducedFactor = 0; % plot toggle

%% Task 01
AR = [4, 6, 8, 10];
b = 10; % [ft]
taper_ratio = linspace(0, 1, 100);
N = 50; % same number of expansions from Anderson

% Initialize values
e = zeros(length(AR), length(taper_ratio));
c_L = zeros(length(AR), length(taper_ratio));
c_Di = zeros(length(AR), length(taper_ratio));

a0_t = 2 * pi; 
a0_r = 2 * pi; 
aero_t = -3; 
aero_r = -3; 
geo_t = 0; 
geo_r = 0;

for i = 1:length(AR)
    
    for j = 1:length(taper_ratio)
        
        taper = taper_ratio(j);
        c_r = (2*b) / (AR(i) * (1 + taper));
        c_t = taper * c_r;
        [e(i,j), c_L(i,j), c_Di(i,j)] = PLLT(b, a0_t, a0_r, c_t, c_r, aero_t, aero_r, geo_t, geo_r, N);

    end
end

if plotInducedFactor
    figure();
    hold on;
    for i = 1:length(AR)
        plot(taper_ratio, (1 ./ e(i, :)) - 1, "LineWidth", 1.5);
    end
    hold off;
    grid on;
    title("Induced Drag Factor as a Function of Taper Ratio", "FontSize", 16);
    ylabel("Induced Drag Factor");
    xlabel("Taper Ratio"); 
    legend("AR 4", "AR 6", "AR 8", "AR 10", "Fontsize", 14);
end

%% Part 03
c_2412=5+(4/12); % chord length (ft)
N=21; % number of panels

%NACA 2412 ----------------
param_2412 = struct("m", 0.02, "p", 0.4, "t", 0.12 * c_2412);
[x_2412, y_2412] = NACA_airfoils(param_2412.m, param_2412.p, param_2412.t, c_2412, N); 

c_0012=3+(8.5/12); % chord length (ft)

%NACA 0012 ----------------
param_0012 = struct("m", 0, "p", 0, "t", 0.12 * c_0012);
[x_0012, y_0012] = NACA_airfoils(param_0012.m, param_0012.p, param_0012.t, c_0012, N); 
alpha = linspace(-10, 10, 200);

[flipped_x_2412, flipped_y_2412] = flipPositions(x_2412, y_2412);
[flipped_x_0012, flipped_y_0012] = flipPositions(x_0012, y_0012);

for i = 1:length(alpha)     
    cl_2412(i) = Vortex_Panel(flipped_x_2412, flipped_y_2412, 1, alpha(i));
    cl_0012(i) = Vortex_Panel(flipped_x_0012,flipped_y_0012, 1, alpha(i));
end

zero_lift_aoa_2412 = interp1(cl_2412, alpha, 0, "linear"); % find alpha for cl = 0 -> L = 0
zero_lift_aoa_0012 = interp1(cl_0012, alpha, 0, "linear"); % find alpha for cl = 0 -> L = 0

lift_slope_2412  = calculateLiftSlope(alpha, cl_2412);
lift_slope_0012 = calculateLiftSlope(alpha, cl_0012);

b_Cessna140 = 33+(4/12);
a0_t_Cessna140 = lift_slope_0012*180/pi; 
a0_r_Cessna140 = lift_slope_2412*180/pi; 
aero_t_Cessna140 = zero_lift_aoa_0012; 
aero_r_Cessna140 = zero_lift_aoa_2412; 
geo_t_Cessna140 = 4; 
geo_r_Cessna140 = 5;
c_t_Cessna140 = c_0012;
c_r_Cessna140 = c_2412;

cessna_N = linspace(1, 300, 300);

for i = 1:length(cessna_N)
    
    odd_term(i) = 2 * i - 1;
    [e_Cessna140(i), c_L_Cessna140(i), c_Di_Cessna140(i)] = PLLT(b_Cessna140, a0_t_Cessna140, a0_r_Cessna140, c_t_Cessna140, c_r_Cessna140, aero_t_Cessna140, aero_r_Cessna140, geo_t_Cessna140, geo_r_Cessna140, i);
   
end

<<<<<<< Updated upstream

alpha=linspace(-16,16,14);
for i=1:length(alpha)
[e_alpha(i), c_L_alpha(i), c_Di_alpha(i)] = PLLT(b_Cessna140, a0_t_Cessna140, a0_r_Cessna140, c_t_Cessna140, c_r_Cessna140, aero_t_Cessna140, aero_r_Cessna140, alpha(i), alpha(i)+1, N);
end

=======
c_L_Cessna140_actual = c_L_Cessna140(length(c_L_Cessna140));
c_Di_Cessna140_actual = c_Di_Cessna140(length(c_Di_Cessna140));

for i = 1:length(cessna_N)
    c_L_Cessna140_error(i) = 100*((c_L_Cessna140(i) - c_L_Cessna140_actual)/c_L_Cessna140_actual);
    c_Di_Cessna140_error(i) = 100*(c_Di_Cessna140(i) - c_Di_Cessna140_actual)/c_Di_Cessna140_actual;
    
    % c_L error
    if c_L_Cessna140_error(i) > 10
        c_L_Cessna140_error_10 = c_L_Cessna140(i+1);
        c_L_error_10_percent = i+1;
    elseif c_L_Cessna140_error(i) > 1
        c_L_Cessna140_error_1 = c_L_Cessna140(i+1);
        c_L_error_1_percent = i+1;
    elseif c_L_Cessna140_error(i) > 0.1
        c_L_Cessna140_error_0_1 = c_L_Cessna140(i+1);
        c_L_error_0_1_percent = i+1;
    end

    % c_Di error
    if c_Di_Cessna140_error(i) > 10
        c_Di_Cessna140_error_10 = c_Di_Cessna140(i+1);
        c_Di_error_10_percent = i+1;
    elseif c_Di_Cessna140_error(i) > 1
        c_Di_Cessna140_error_1 = c_Di_Cessna140(i+1);
        c_Di_error_1_percent = i+1;
    elseif c_Di_Cessna140_error(i) > 0.1
        c_Di_Cessna140_error_0_1 = c_Di_Cessna140(i+1);
        c_Di_error_0_1_percent = i+1;
    end
end

figure()
plot(cessna_N,c_L_Cessna140_error)
hold on
yline(10)
yline(1)
yline(0.1)










>>>>>>> Stashed changes
