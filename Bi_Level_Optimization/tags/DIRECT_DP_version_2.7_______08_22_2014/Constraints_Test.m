% clear all
% close all
% clc
% tic 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%----------------------------Load All Data--------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

cd('Components');
%                              ~~ Engine ~~
Engine_2rz_0410;   % Set all optimal engine speeds
% Engine_102_kW;
% Engine_41_kW;

%                              ~~ Motor ~~
% Motor_int;
% Motor_75_kW;
Motor_30_kW;

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
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%----------------Update the Design Variables------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
dvar.FD = 3.65;
dvar.G = 1.7;
dvar.fc_trq_scale =  0.7;
dvar.mc_trq_scale = 0.8;
mc_max_pwr_kW =  dvar.mc_trq_scale*vinf.mc_max_pwr_kW;
% dvar.module_number = ceil(4*mc_max_pwr_kW*1000*Rint_size/(Voc_size^2));
dvar.module_number = 15;  % Fixed (for now) - also NEED to change in both obj and con!!!!!
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%-----------Manipulate Data Based of Scaling Factors----------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
Manipulate_Data_Structure; % Need to recalcualte the Tw for the ne vehicle mass

% TEMP___
cd('Drive_Cycle')
Manipulate_Drive_Cycle;
cd ..
break;
% TEMP___

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%-------------------------Acceleration Tests ------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Initialize Stuff
V_sim_full = [];
We_sim_full = [];
Wm_sim_full = [];
x_sim_full = [];
acc_sim_full = [];
Teng_sim_full = [];
Tm_sim_full = [];
time_sim_full = [];

cd('Initial Component Sizing')

index(1) = [1];
n = 1;
V_0 = 0;
V_f = 60;
dt_2 = 12;
Acc_Final_new = 100;  % Does not matter
TYPE = 1; % Velocity req.
[ pass_acc_test(n), Sim_Variables ] = Acceleration_Test(V_0,V_f, Acc_Final_new, dt_2, param, vinf, dvar, TYPE);
index(n+1) = [index(n) + length(Sim_Variables(2:end,1))];

% Save Variables
x_sim_full = [ x_sim_full;Sim_Variables(:,1)];
V_sim_full = [V_sim_full; Sim_Variables(:,2)]; % MPH
acc_sim_full = [ acc_sim_full; Sim_Variables(:,3)];
Teng_sim_full = [Teng_sim_full; Sim_Variables(:,4)];
Tm_sim_full = [Tm_sim_full; Sim_Variables(:,5)];
We_sim_full = [We_sim_full; Sim_Variables(:,6)];   % RPM
Wm_sim_full = [ Wm_sim_full; Sim_Variables(:,7)];   % RPM
time_sim_full = [time_sim_full; Sim_Variables(:,8)];

dt_2 = 0.0002;
load V_0_new;
load V_f_new;
load Acc_Final_new
TYPE = 0; % Acceleration req.
for i = 1:length(V_0_new)
n = n + 1;
[ pass_acc_test(n), Sim_Variables ] = Acceleration_Test(V_0_new(i),V_f_new(i), Acc_Final_new(i),dt_2, param, vinf, dvar, TYPE);
index(n+1) = [index(n) + length(Sim_Variables(2:end,1))];

% Save Variables
x_sim_full = [ x_sim_full;Sim_Variables(:,1)];
V_sim_full = [V_sim_full; Sim_Variables(:,2)]; % MPH
acc_sim_full = [ acc_sim_full; Sim_Variables(:,3)];
Teng_sim_full = [Teng_sim_full; Sim_Variables(:,4)];
Tm_sim_full = [Tm_sim_full; Sim_Variables(:,5)];
We_sim_full = [We_sim_full; Sim_Variables(:,6)];   % RPM
Wm_sim_full = [ Wm_sim_full; Sim_Variables(:,7)];   % RPM
time_sim_full = [time_sim_full; Sim_Variables(:,8)];
end

fail_acc_test = ~pass_acc_test;
FAIL_ACCEL_TEST = any(fail_acc_test)

e = 3;
start = index(e)+ e;
stop = index(e+1);
x_sim = x_sim_full(start:stop);
V_sim =  V_sim_full(start:stop);
acc_sim = acc_sim_full(start:stop);
Teng_sim =  Teng_sim_full(start:stop);
Tm_sim = Tm_sim_full(start:stop);
We_sim = We_sim_full(start:stop);
Wm_sim =  Wm_sim_full(start:stop);
time_sim = time_sim_full(start:stop);

% Plot the results
Kinematics_Plot;
Kinetics_Plot;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%-----------------------------Grade Test----------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%--------------------------Set Requirements--------------------------------
Motor_ON = 1;

% Test 1
r = 2;
V_test(r) = 80*param.mph_mps;
alpha_test(r) = 0*pi/180;

% Test 2
r = 1;
V_test(r) = 55*param.mph_mps;
alpha_test(r) = 5*pi/180;

[Sim_Grade, FAIL_GRADE_TEST] = Grade_Test( param, vinf, dvar, alpha_test, V_test, Motor_ON );
RR = length(alpha_test);
FAIL_GRADE_TEST
Grade_Plot;
cd ..