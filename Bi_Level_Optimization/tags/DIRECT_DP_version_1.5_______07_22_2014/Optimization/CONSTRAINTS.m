function [con, con_e] = CONSTRAINTS(x,varargin) 

con = evalin('base','con');

con_e=0;

return
