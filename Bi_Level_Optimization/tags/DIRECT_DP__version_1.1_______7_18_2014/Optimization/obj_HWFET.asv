function obj=obj_HWFET(x,varargin)

% initialize
% error=0;
obj=0;

% Reloading Cycle and Other Parameters that may have been changed by grade
% and accel test
% temp_cyc_mph = evalin('base', 'cyc_mph_original');
% assignin('base', 'cyc_mph', temp_cyc_mph);
% temp_time_final = evalin('base', 'time_final_original');
% assignin('base', 'time_final', temp_time_final);
% assignin('base', 'theta', 0);

% Update the Design Variables
assignin('base','FD',x(1));
% assignin('base','mc_trq_scale',x(2));
% assignin('base','ess_module_num',x(3));
% assignin('base','fd_ratio',x(4));
% assignin('base','ess_cap_scale',x(5));
% assignin('base','M_ratio',x(6));
% assignin('base','E_ratio',x(7));
% 
% fc_trq_scale = x(1);
% mc_trq_scale = x(2);
% fd_ratio = x(4);
% M_ratio = x(6);
% E_ratio = x(7);
% 
% Create_Motor_Shift_Map;
% Create_Engine_Shift_Map;
% assignin('base','M_up',M_up);
% assignin('base','M_down',M_down);
% assignin('base','M_up_mot',M_up_mot);
% assignin('base','M_down_mot',M_down_mot);

% gear_ratio = evalin('base','gear_ratio');
% Constant = evalin('base', 'Constant');
% V_limit = Constant./gear_ratio./fd_ratio;
% assignin('base', 'V_limit', V_limit);
% V_max_limit = ceil(V_limit(end))+5;
% assignin('base', 'V_max_limit', V_max_limit);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Mass Scaling
% ess_mass=inline('(f(1)*ess_module_num+f(2))*(f(3)*ess_cap_scale+f(4))*(ess_module_mass)','f','ess_module_num','ess_cap_scale','ess_module_mass');
% ess_mass_coef=[1 0 1 0]; % use as f
% ess_module_mass = evalin('base','ess_module_mass');
% 
% init_fc_trq_scale = evalin('base', 'init_fc_trq_scale');
% init_mc_trq_scale = evalin('base', 'init_mc_trq_scale');
% 
% mass_engine = 215.0776*fc_trq_scale/init_fc_trq_scale;
% mass_motor = 60*mc_trq_scale/init_mc_trq_scale; 
% mass_bat = ess_mass(ess_mass_coef,x(3),x(5),ess_module_mass);
% 
% comp_mass = mass_engine + mass_motor + mass_bat;
% base_mass = evalin('base','base_mass');
% passenger_mass = evalin('base','passenger_mass');
% assignin('base', 'veh_glider_mass', comp_mass+base_mass+passenger_mass);
% 
% Je_engine = (1.9/(0.224809*12*3.28))*fc_trq_scale;
% assignin('base', 'Je_engine', Je_engine);
% 
% m_inertia = 0.0507*mc_trq_scale;
% assignin('base', 'm_inertia', m_inertia);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get_Params;
FD = x(1);
G = x(2);

% Run the model in the base workspace
heidelberg_DIRECT;

% 
% FAIL =  evalin('base', 'FAIL'); % Do I need this?

if FAIL == 1; % It failed
    FAIL_DP = 100;
else
    FAIL_DP = 0.5;
end

    
%     if fail == 0;

% Huck added delta SOC
%         delta_soc =  evalin('base', 'delta_soc');
%         delta_trace =  evalin('base', 'delta_trace');;
% 
% MPG =  evalin('base', 'MPG');

%         CO =  evalin('base', 'CO');
%         NOx =  evalin('base', 'NOx');
%         HC =  evalin('base', 'HC');
%         PM =  evalin('base', 'PM');
%
%         mpg_nominal = evalin('base', 'mpg_nominal');
%         CO_nominal = evalin('base', 'CO_nominal');
%         HC_nominal = evalin('base', 'HC_nominal');
%         PM_nominal = evalin('base', 'PM_nominal');
%         NOx_nominal = evalin('base', 'NOx_nominal');

%         SC_bat = evalin('base', 'SC_bat');
%         Batt_eff_discharge = evalin('base', 'Batt_eff_discharge');
%         Batt_eff_charge = evalin('base', 'Batt_eff_charge');
%         Eff_elec_mot = evalin('base', 'Eff_elec_mot');

obj = MPG; % UPDATE WITH EMMISSIONS/PERFORMANCE

%     else  % Fail it!!!
%         delta_soc = 100; delta_trace = 100; MPG= 0; CO = 0; HC = 0; NOx = 0; PM = 0; SC_bat = 0; Batt_eff_discharge = 0; Batt_eff_charge = 0; Eff_elec_mot = 0;
%     end

assignin('base','con',FAIL_DP);

return
