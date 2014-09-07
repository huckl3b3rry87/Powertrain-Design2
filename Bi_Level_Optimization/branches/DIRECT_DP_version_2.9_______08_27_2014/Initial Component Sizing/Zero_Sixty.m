%% Initialize Stuff
V_sim_full = [];
We_sim_full = [];
Wm_sim_full = [];
x_sim_full = [];
acc_sim_full = [];
Teng_sim_full = [];
Tm_sim_full = [];
time_sim_full = [];

dt_2 = 12;
index(1) = [1];

V_0 = 0;
V_f = 60*param.mph_mps;

[ PASS, Sim_Variables ] = Acceleration_Test(V_0,V_f, dt_2, param, vinf, dvar);
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