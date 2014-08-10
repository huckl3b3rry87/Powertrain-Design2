% main.m
% A sample program to solve Constant
% Power Performance (CPP)
% of vehicles in longitudinal motion
clc , close all , clear all
global P c F0 m % Shares these information with
% function const_pow
% Vehicle information:
m=1000; % Vehicle mass (kg)
c=0.4; % Aerodynamics overall coefficient
fr=0.02; % Rolling resistance coefficient
P=60000; % Input power (watts)
theta=0; % Angle of slope (rad)

% Determine the total constant resistive force:
F0=(fr*cos(theta)+sin(theta))*m*9.81;
% Initial condition:
% Initial velocity cannot be 'zero' due to 'division
% by zero' in P/v term
v0=eps;
% define differentiation interval t0-tf:
t0=0; tf=10;
% Now invoke ode45:
[t,v]=ode45(@const_pow, [t0 tf], v0);
% Calls ‘const_pow’ function
% Plot the variation of velocity versus time
plot(t,v)
xlabel('Time (s)')
ylabel('Velocity (m/s)')
grid

