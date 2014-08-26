clear all
clc
close all

load V_0;
load V_f;

% Later this will be done by DIRECT

dvar.fc_trq_scale = 1.0982;
dvar.mc_trq_scale = 1.5;
dvar.FD = 4.4800;
dvar.G = 1.4000;
dvar.module_number = 50;

% Enter proper vehicle parameters:
cd('Components');

% Eninge
Engine_41_kW;
% Engine_102_kW;

% Motor
Motor_75_kW;

module_number = 50;  % Fix later
Battery_ADVISOR;

% Calcualtes New Vehicle Mass
Vehicle_Parameters_4_HI_AV; % mass of motor will change in here..
cd ..

data; % Put all the data into structures

% Later the data will be passed to this function

% Initialize Stuff
V_sim_full = [];
We_sim_full = [];
Wm_sim_full = [];
x_sim_full = [];
acc_sim_full = [];
Teng_sim_full = [];
Tm_sim_full = [];
time_sim_full = [];

dt_2 = 0.1;
index(1) = [1];
for i = 1:length(V_0)
    [ PASS, Sim_Variables ] = Acceleration_Test(V_0(i),V_f(i), dt_2, param, vinf, dvar);
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
%%
% e = 6;

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

