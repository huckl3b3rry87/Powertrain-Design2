clear all
clc
close all
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%------------------Load the Component Data--------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
cd('Components');

% Engine
% Engine_41_kW;
% Engine_102_kW;
Engine_2rz_0410;

% Motor
% Motor_75_kW;
Motor_30_kW;

Battery_ADVISOR;

% Vehicle
Vehicle_Parameters_4_HI_AV;
cd ..

data;                              %  Put all the data into structures
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%-------------------------Design Variables--------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
dvar.fc_trq_scale = 0.8700;
dvar.mc_trq_scale = 1;
dvar.FD = 3.65;
dvar.G = 1.6;  % Also change this so that you can get to the maximum speed!

mc_max_pwr_kW =  dvar.mc_trq_scale*vinf.mc_max_pwr_kW;
dvar.module_number = ceil(4*mc_max_pwr_kW*1000*Rint_size/(Voc_size^2));
Manipulate_Data_Structure;     % Update the data based off of the new variables

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%---Determine Acceleration Req. Based off Drive Cycles--------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
clear V_0 V_f
S = 10;                    % Number of points to calcualte acceleration at
dt_2 = 0.0002;
graph = 0;
[ V_0, V_f, Acc_Final ] = Accel_Req_2( S, dt_2, graph );  % V_Final is in MPH
% [ V_0, V_f, Acc_Final ] = Accel_Req( S, dt_2, graph );  % Includes US06
% Initialize Stuff
V_sim_full = [];
We_sim_full = [];
Wm_sim_full = [];
x_sim_full = [];
acc_sim_full = [];
Teng_sim_full = [];
Tm_sim_full = [];
time_sim_full = [];

dvar.mc_trq_scale = 1;
mc_max_pwr_kW =  dvar.mc_trq_scale*vinf.mc_max_pwr_kW;
dvar.module_number = ceil(4*mc_max_pwr_kW*1000*Rint_size/(Voc_size^2));
Manipulate_Data_Structure;

index(1) = [1];
for i = 1:(length(V_0)+1)
    
    if i == (length(V_0)+1)
        TYPE = 1; % Velocity req.
        V_0_n = 0;
        V_f_n = 60;
        dt_2 = 12;
        Acc_Final = 100;  % Does not matter
        [ PASS, Sim_Variables ] = Acceleration_Test(V_0_n,V_f_n, Acc_Final, dt_2, param, vinf, dvar, TYPE);
    else
        TYPE = 0; % Acceleration req.
        dt_2 = 0.0002;
        [ PASS, Sim_Variables ] = Acceleration_Test(V_0(i),V_f(i), Acc_Final(i), dt_2, param, vinf, dvar, TYPE);
        
    end
    Pass_Total(i) = PASS;
    index(i+1) = [index(i) + length(Sim_Variables(2:end,1))];
    
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

fprintf('Did the vehicle pass??? \n\n')
I = [];
if any(Pass_Total == 0)
    fprintf('Nope!! \n')
    I = find(Pass_Total == 0)
    fprintf('\n\n')
else
    fprintf('Yep!! \n\n')
end

if ~isempty(I)
    e = I(1);
else
    e = 1;
end
%
e = (length(V_0)+1);
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





