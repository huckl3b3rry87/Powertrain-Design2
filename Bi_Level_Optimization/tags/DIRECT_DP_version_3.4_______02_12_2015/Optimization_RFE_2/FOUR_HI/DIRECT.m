clear all
close all
clc
tic 

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%----------------------------Load All Data--------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

cd('Components');
%                              ~~ Engine ~~
% Engine_73_kW;   
% Engine_102_kW;
% Engine_41_kW;
Engine_50_kW;

%                              ~~ Motor ~~
% Motor_int;
% Motor_75_kW;
% Motor_30_kW;
Motor_25_kW;
% Motor_49_kW;
% Motor_16_kW;

%                             ~~ Battery ~~
% Battery_int;  % No variation with the number of modules in this battery!!
Battery_ADVISOR;

%                              ~~ Vehicle ~~

% Vehicle_Parameters_4_HI_AV;
Vehicle_Parameters_4_HI;

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
e = 1;  % If you change engine and motor get ride of e on fc_..
dvar.FD = 5.655*e;
dvar.G = 1.1*e;
dvar.fc_trq_scale =  1.015*e + vinf.Paux/(vinf.fc_max_pwr*1000);
dvar.mc_trq_scale = 1.2*e;
mc_max_pwr_kW =  dvar.mc_trq_scale*vinf.mc_max_pwr_kW;
% dvar.module_number = ceil(4*mc_max_pwr_kW*1000*Rint_size/(Voc_size^2));
dvar.module_number = 20*e;  % Fixed (for now) - also NEED to change in both obj and con!!!!!

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%---------------------Select Drive Cycle----------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%                              ~~ Standard ~~
% cyc_name = 'HWFET_AA';
cyc_name = 'HWFET';
% cyc_name = 'UDDS';
% cyc_name = 'US06';
% cyc_name = 'SHORT_CYC_HWFET';
% cyc_name = 'RAMP';
% cyc_name = 'LA92';
% cyc_name = 'CONST_65';
% cyc_name = 'CONST_45';
% cyc_name = 'COMMUTER';


%                              ~~ AV~~

% cyc_name = 'US06_AV';
% cyc_name = 'HWFET_AV';

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%---------------------Run Optimization------------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
L_eng = 0.2;
H_eng = 1.8;
L_g = 0.1;
H_g = 1.9;

% Identify the Design Variables and their ranges    
dv_names={ 'FD', 'G','fc_trq_scale','mc_trq_scale'};
x_L=[    L_g*dvar.FD, L_g*dvar.G, L_eng*dvar.fc_trq_scale, L_eng*dvar.mc_trq_scale]';
x_U=[    H_g*dvar.FD, H_g*dvar.G, H_eng*dvar.fc_trq_scale, H_eng*dvar.mc_trq_scale]';
   
con_names={'FAIL', 'delta_SOC','FAIL_ACCEL_TEST','FAIL_GRADE_TEST'}; % They Also Did not Do Grade and Acceleration Tests!
c_L=[          0;      -0.01;         0;                0];
c_U=[          0;       0.01;         0;                0];

% Define the Problem
PriLev = 2;                      % 1 is no graph 2 shows a graph
MaxEval = 500;
MaxIter = inf;
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

% start the optimization                                                                          %    1          2          3          4            5             6          7         8         
DP_optimization = gclSolve(p_f, p_c, x_L, x_U, A, b_L, b_U, c_L, c_U, I, GLOBAL, PriLev, plot_info, dv_names, resp_names, con_names, param_names, param_data, vinf_names, vinf_data, cyc_name );

% save the results
eval(['save(''',prev_results_filename,''',','''DP_optimization'');']) 

PlotOptimResults(DP_optimization.GLOBAL.f_min_hist, plot_info)

% end timer
toc