% clear all
% clc
% close all

load V_0;
load V_f;

% Later this will be done by DIRECT

fc_trq_scale = 1.0982;
mc_trq_scale = 1.2656;
FD = 4.4800;
G = 1.4000;

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

dt_2 = 0.1;

for i = 1:length(V_0)
    
% [ PASS, Sim_Variables ] = Acceleration_Test(V_0(i),V_f(i), dt_2);

evalin('base','[ PASS, Sim_Variables ] = Acceleration_Test(V_0(i),V_f(i), dt_2)')

Pass_Total(i) = PASS;

end

fprintf('Did the vehicle pass?')

any(Pass_Total == 0)
