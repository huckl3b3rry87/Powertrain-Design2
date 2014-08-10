% Function called in main program ‘main.m’
function f=const_pow(t,v)
global P c F0 m % Shares these information
% with main program
% define differential equation to be solved
% (f=dv/dt)
f=(P/v-F0-c*v^2)/m;