clear all
close all
clc
tic 

FD = 3.4; % Final Drive Ratio
G = 3;
% %% Pass the vehicle
% heidelberg_DIRECT_SETUP;

% drivingcycle = 2; %1:EPA HWFET cycle %2:EPA Highway cycle
% SETUP;
% load Project_eng_params;
% Get_Params;
% 
% temp_time_final = evalin('base', 'time_final_original');
% temp_cyc_mph = evalin('base', 'cyc_mph_original');
% actual_spd = Results.get('actual_spd');
% samp_time = evalin('base', 'samp_time');
% delta_spd = zeros(length(temp_time_final));
% for i = 1:temp_time_final
%     delta_spd(i) = temp_cyc_mph(i,2) - actual_spd(i/samp_time -9); % Assume it is at the end
% end
% 
% delta_trace_nominal = max(abs(delta_spd));
% delta_soc_nominal = delta_soc;
% mpg_nominal = Results.get('mpg');
% mpg_nominal = mpg_nominal(length(mpg_nominal))
% CO_nominal = Results.get('CO');
% CO_nominal = CO_nominal(length(CO_nominal));
% HC_nominal = Results.get('HC');
% HC_nominal = HC_nominal(length(HC_nominal));
% NOx_nominal = Results.get('NOx');
% NOx_nominal = NOx_nominal(length(NOx_nominal));
% PM_nominal = Results.get('PM');
% PM_nominal = PM_nominal(length(PM_nominal));

%%%%%%%%%%% RUN OPTIMIZATION %%%%%%%%%%%%%%%%%%%%%%

% Identify the Design Variables and their ranges    
dv_names={ 'FD', 'G'};
x_L=[    0.5*FD, 0.5*G]';
x_U=[    1.5*FD, 0.5*G]';
   
con_names={'FAIL_DP'}; % They Also Did not Do Grade and Acceleration Tests!
c_L=[        0;     ];
c_U=[        1;     ];

% Define the Problem
PriLev = 2;                      % 1 is no graph 2 shows a graph
MaxEval = 50;
MaxIter = 49;
GLOBAL.epsilon = 1e-4;

% Define the Objective function Name for the GRAPH
resp_names={'DP_HWFET'};

p_f='obj_HWFET';
p_c='con_HWFET';

A=[];
b_L=[];
b_U=[];
I=[];
cont_bool=0;    % Continue from a Previous Run??
prev_results_filename='DP_HWFET';

if cont_bool==1
   eval(['load(''',prev_results_filename,''')']) 
   GLOBAL = DP_HWFET_optimization.GLOBAL;
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

% start the optimization
DP_HWFET_optimization = gclSolve(p_f, p_c, x_L, x_U, A, b_L, b_U, c_L, c_U, I, GLOBAL, PriLev, plot_info, dv_names, resp_names, con_names);

% save the results
eval(['save(''',prev_results_filename,''',','''DP_HWFET_optimization'');']) 

PlotOptimResults(DP_HWFET_optimization.GLOBAL.f_min_hist, plot_info)

% end timer
toc