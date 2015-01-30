function [con, con_e] = CONSTRAINTS(x,varargin) 

% FAIL and delta_SOC come in  here
con = evalin('base','con');
offset=length(con);

param =  cell2struct(varargin{5}, varargin{4},1);
vinf =  cell2struct(varargin{7}, varargin{6},1);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%----------------Update the Design Variables------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
dvar.FD = x(1);
dvar.G = x(2);
dvar.fc_trq_scale = x(3);
dvar.mc_trq_scale = x(4);  
dvar.module_number = 12;  % Fixed (for now) - should be passing this..
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%-----------Manipulate Data Based of Scaling Factors----------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
Manipulate_Data_Structure; % Need to recalcualte the Tw for the ne vehicle mass

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%-------------------------Acceleration Tests ------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
cd('Initial Component Sizing')

n = 1;
V_0 = 0;
V_f = 60;
dt_2 = 12;
Acc_Final_new = 100;  % Does not matter
TYPE = 1; % Velocity req.
[ pass_acc_test(n), Sim_Variables, time_sixty, Acc_Test] = Acceleration_Test(V_0,V_f, Acc_Final_new, dt_2, param, vinf, dvar, TYPE);

dt_2 = 0.0002;
load V_0;
load V_f;
load Acc_Final;
TYPE = 0; % Acceleration req.
for i = 1:length(V_0_new)
    n = n + 1;
    [ pass_acc_test(n), Sim_Variables, time_sixty, Acc_Test] = Acceleration_Test(V_0_new(i),V_f_new(i), Acc_Final_new(i),dt_2, param, vinf, dvar, TYPE);
end

fail_acc_test = ~pass_acc_test;
FAIL_ACCEL_TEST = any(fail_acc_test);

if ~isempty(FAIL_ACCEL_TEST)  % &~isempty(z0_60)&~isempty(z0_85)
    con(offset+1,1)= FAIL_ACCEL_TEST;
    %    con(offset+2,1)=time40_60;
    %    con(offset+3,1)=time0_85;
else
    %    con(offset+1:offset+3,1)=100;
    con(offset+1,1)= 0;
end

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

if ~isempty(FAIL_GRADE_TEST)
    con(offset+2,1)= FAIL_GRADE_TEST;
else
    con(offset+2,1) = 0;
end

con_e=0;
cd .. 
return
