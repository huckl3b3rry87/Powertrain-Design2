function obj = OBJECTIVE(x,varargin)

param =  cell2struct(varargin{5}, varargin{4},1);
vinf =  cell2struct(varargin{7}, varargin{6},1);
cyc_name = varargin{8};
RUN_TYPE = cell2struct(varargin{10},varargin{9},1);
weight = cell2struct(varargin{12},varargin{11},1);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%----------------Update the Design Variables------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
dvar.FD = x(1);
dvar.G = x(2);
dvar.fc_trq_scale = x(3);
dvar.mc_trq_scale = x(4);  
dvar.module_number = 38;  % Fixed (for now) - should be passing this..
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%-----------Manipulate Data Based of Scaling Factors----------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
Manipulate_Data_Structure; % Need to recalcualte the Tw for the ne vehicle mass

[cyc_data] = Drive_Cycle(param, vinf, cyc_name );

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%---------------------Run DP with new Data--------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%    

[FAIL, MPG, emission, delta_SOC, sim] = Dynamic_Programming_func(param, vinf, dvar, cyc_data, RUN_TYPE, weight);

obj = -MPG + weight.NOx*emission.NOx + weight.CO*emission.CO + weight.HC*emission.HC; % UPDATE WITH EMMISSIONS/PERFORMANCE

assignin('base','con',[FAIL.final; delta_SOC; MPG; emission.NOx; emission.CO; emission.HC]);

return
