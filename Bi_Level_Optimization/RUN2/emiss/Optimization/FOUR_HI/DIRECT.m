clear all
close all
clc
tic

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%----------------------------Load All Data--------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

cd('Components');
%                              ~~ Engine ~~
% Engine_30_kW;
Engine_41_kW;
% Engine_41_kW_smooth;
% Engine_50_kW;
% Engine_73_kW;
% Engine_95_kW;
% Engine_102_kW;
% Engine_186_kW;
% Engine_224_kW;

%                              ~~ Motor ~~
% Motor_int;
% Motor_75_kW;
% Motor_30_kW;
Motor_49_kW;
% Motor_10_kW;
% Motor_8_kW;
% Motor_16_kW;

%                             ~~ Battery ~~
% Battery_int;  % No variation with the number of modules in this battery!!
Battery_ADVISOR;

%                              ~~ Vehicle ~~
% Vehicle_Parameters_4_HI_AV;
Vehicle_Parameters_4_HI_small;
% Vehicle_Parameters_8_HI_AV;
% Vehicle_Parameters_8_HI;

% Low Speed
% Vehicle_Parameters_4_low_AV;
% Vehicle_Parameters_4_low;
% Vehicle_Parameters_1_low_AV;
% Vehicle_Parameters_1_low;

cd ..
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%-------------Put all the data into structures and cells------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
data;     
param_names = fieldnames(param);
param_data = struct2cell(param);
vinf_names = fieldnames(vinf);
vinf_data = struct2cell(vinf);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%---------------------Update the Design Variables-------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~d%
e = 1;
dvar.FD = 5.57*e;
dvar.G = 1.3*e;
dvar.fc_trq_scale = 0.895;
dvar.mc_trq_scale = 1.2*e;
mc_max_pwr_kW =  dvar.mc_trq_scale*vinf.mc_max_pwr_kW;
% dvar.module_number = ceil(4*mc_max_pwr_kW*1000*Rint_size/(Voc_size^2));
dvar.module_number = 38;% Fixed (for now) - also NEED to change in both obj and con!!!!!
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%---------------------Update the Data-------------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
Manipulate_Data_Structure;  % May not have to do this here
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%---------------------Select Drive Cycle----------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%%                              ~~ Standard ~~

% cyc_name = 'HWFET';
cyc_name = 'UDDS';
% cyc_name = 'US06';
% cyc_name = 'SHORT_CYC_HWFET';
% cyc_name = 'RAMP';
% cyc_name = 'LA92';
% cyc_name = 'CONST_65';
% cyc_name = 'CONST_45';
% cyc_name = 'COMMUTER';

% City
% cyc_name = 'INDIA_URBAN';
% cyc_name = 'MANHATTAN';
% cyc_name = 'Nuremberg';
% cyc_name = 'NYCC';
% cyc_name = 'AA_final';
%                              ~~ AV~~

% cyc_name = 'US06_AV';
% cyc_name = 'HWFET_AV';
% cyc_name = 'AA_final_AV';

[cyc_data] = Drive_Cycle(param, vinf, cyc_name );

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%-----------------Define the Run Type-------------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
RUN_TYPE.sim = 1;  % RUN_TYPE = 1 - for DIRECT     &    RUN_TYPE = 0 - for DP only
RUN_TYPE.emiss = 1; % RUN_TYPE.emiss = 1 - has emissions  &   RUN_TYPE.emiss = 0 - NO emissions
RUN_TYPE.plot = 0;  % RUN_TYPE.plot = 1 - plots on  &   RUN_TYPE.plot = 0 - plots off
RUN_TYPE.soc_size = 0.005;
RUN_names = fieldnames(RUN_TYPE);
RUN_data = struct2cell(RUN_TYPE);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%-----------------Weighing Parameters for DP------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
weight.fuel = 1*1.4776/1.4776;  % These are for a specific engine, we need to change this!
if RUN_TYPE.emiss == 1
    weight.NOx = 0.1*1.4776/0.0560;
    weight.CO = 0.6*1.4776/0.6835;
    weight.HC = 0.3*1.4776/0.0177;
end
weight.shift = 0.3;
weight.engine_event = 10; % Is then multiplied by fc_trq_scale
weight_names = fieldnames(weight);
weight_data = struct2cell(weight);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%---------------------Run Optimization------------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%


% Identify the Design Variables and their ranges    
dv_names={ 'FD', 'G','fc_trq_scale','mc_trq_scale'};
x_L=[    0.5*dvar.FD, 0.5*dvar.G, 0.5*dvar.fc_trq_scale, 0.5*dvar.mc_trq_scale]';
x_U=[    1.5*dvar.FD, 1.5*dvar.G, 1.5*dvar.fc_trq_scale, 1.5*dvar.mc_trq_scale]';
   
con_names={'FAIL',         'delta_SOC',       'MPG',  'NOx',    'CO',  'HC',  'FAIL_ACCEL_TEST', 'FAIL_GRADE_TEST', }; % They Also Did not Do Grade and Acceleration Tests!
c_L=[          0;      -RUN_TYPE.soc_size;      0;      0;        0;     0;           0;                 0];
c_U=[          0;       RUN_TYPE.soc_size;     1000;  1000;     1000;   1000;         0;                 0];

% Define the Problem
PriLev = 2;                      % 1 is no graph 2 shows a graph
MaxEval = 500;
MaxIter = 499;
GLOBAL.epsilon = 1e-4;

% Define the Objective function Name for the GRAPH
resp_names={'DP'};

p_f='OBJECTIVE';
p_c='CONSTRAINTS';

A=[];
b_L=[];
b_U=[];
I=[];
cont_bool=0;    % Continue from a Previous Run??
prev_results_filename='DP';

if cont_bool==1
   eval(['load(''',prev_results_filename,''')']) 
   GLOBAL = DP_optimization.GLOBAL;
   GLOBAL.MaxEval = MaxEval;
   GLOBAL.MaxIter = MaxIter;
else
   GLOBAL.MaxEval = MaxEval;
   GLOBAL.MaxIter = MaxIter;
end

plot_info.var_label=dv_names;
plot_info.var_ub=num2cell(x_U);
plot_info.var_lb=num2cell(x_L);
plot_info.con_label=con_names;
plot_info.con_ub=num2cell(c_U);
plot_info.con_lb=num2cell(c_L);
plot_info.fun_label=resp_names;

% start the optimization                                                                          %    1          2          3          4            5             6          7         8          9         10           11          12       
DP_optimization = gclSolve(p_f, p_c, x_L, x_U, A, b_L, b_U, c_L, c_U, I, GLOBAL, PriLev, plot_info, dv_names, resp_names, con_names, param_names, param_data, vinf_names, vinf_data, cyc_name, RUN_names, RUN_data, weight_names, weight_data );

% save the results
eval(['save(''',prev_results_filename,''',','''DP_optimization'');']) 

PlotOptimResults(DP_optimization.GLOBAL.f_min_hist, plot_info)

% end timer
toc