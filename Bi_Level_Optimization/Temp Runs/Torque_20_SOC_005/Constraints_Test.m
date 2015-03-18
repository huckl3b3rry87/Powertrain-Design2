clear all
close all
clc
tic 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%----------------------------Load All Data--------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

cd('Components');
%                              ~~ Engine ~~
% Engine_2rz_0410;   % Set all optimal engine speeds
% Engine_102_kW;
% Engine_41_kW;
Engine_50_kW;

%                              ~~ Motor ~~
% Motor_int;
% Motor_75_kW;
% Motor_30_kW;
Motor_25_kW;


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
dvar.FD = 5.655;
dvar.G = 1.1;
dvar.fc_trq_scale =  1;
dvar.mc_trq_scale = 1.25;
mc_max_pwr_kW =  dvar.mc_trq_scale*vinf.mc_max_pwr_kW;
% dvar.module_number = ceil(4*mc_max_pwr_kW*1000*Rint_size/(Voc_size^2));
dvar.module_number = 21;  % Fixed (for now) - also NEED to change in both obj and con!!!!!
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%-----------Manipulate Data Based of Scaling Factors----------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
Manipulate_Data_Structure; % Need to recalcualte the Tw for the ne vehicle mass

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
V_0_i = 0;
V_f_i = 60;
dt_2 = 12;
Acc_Final_new = 100;  % Does not matter
TYPE = 1; % Velocity req.
[ pass_acc_test(n), Sim_Variables, time_sixty, Acc_Test  ] = Acceleration_Test(V_0_i,V_f_i, Acc_Final_new, dt_2, param, vinf, dvar, TYPE);
time_sixty_save = time_sixty;

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
load V_0;
load V_f;
load Acc_Final;
TYPE = 0; % Acceleration req.
for i = 1:length(V_0)
n = n + 1;
[ pass_acc_test(n), Sim_Variables, time_sixty, Acc_Test_ ] = Acceleration_Test(V_0(i),V_f(i), Acc_Final(i),dt_2, param, vinf, dvar, TYPE);
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

% time_sixty_save{i+1} = time_sixty;
i
Acc_Test_save{i} = Acc_Test_

end

fail_acc_test = ~pass_acc_test;
FAIL_ACCEL_TEST = any(fail_acc_test)

e = 1;
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
Motor_ON = 0;

% Test 1
r = 1;
V_test(r) = 80*param.mph_mps;
alpha_test(r) = 0*pi/180;

% Test 2
r = 2;
V_test(r) = 55*param.mph_mps;
alpha_test(r) = 5*pi/180;

[Sim_Grade, FAIL_GRADE_TEST, V_max] = Grade_Test( param, vinf, dvar, alpha_test, V_test, Motor_ON );
RR = length(alpha_test);
FAIL_GRADE_TEST
Grade_Plot;
cd ..