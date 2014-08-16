function [con, con_e] = CONSTRAINTS(x,varargin) 

con = evalin('base','con');

% Acceleration and Grade Tests

% Delta SOC



con_e=0;

return
