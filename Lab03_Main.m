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

%Task 1
% 
% m= % maximum chamber (percent of chord)
% p= % x location of maximum chamber (in thenths of chord)
% t= % maximum thickness (percent of chord)


c=1 % chord length (m)
N=50 % number of panels
n=1000 %number of x coordinates

%% task 2
m = 0; 
p = 0;
t = 12/100;
alpha = 12; % degrees
v_inf = 50;
error = 100;
c_l_actual = 1;
N = 100;

while error > 0.01
    N = N + 1;
    [x_b,y_b] = NACA_airfoils(m,p,t,c,N);
    [CL,CP,CIRC,X,Y] = Vortex_Panel_2(x_b,y_b,v_inf,alpha,flag);
    error = (CL - c_l_actual)/c_l_actual;
end


