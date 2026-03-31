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

m=[0;2]./100; % maximum chamber (percent of chord)
p=[0;4]./10; % x location of maximum chamber (in thenths of chord)
t=[21;21]./100; % maximum thickness (percent of chord)


c=1; % chord length (m)
N=50; % number of panels

x_geo=Geometric_Distribution(c,N);

%NACA 0021

[x_b_0021, y_b_0021] = NACA_airfoils(m(1),p(1),t(1),c,N);



%NACA 2421
[x_b_2421, y_b_2421] = NACA_airfoils(m(2),p(2),t(2),c,N);


%Checking
figure
hold on
plot(x_b_0021,y_b_0021)
plot(x_b_2421,y_b_2421)
hold off