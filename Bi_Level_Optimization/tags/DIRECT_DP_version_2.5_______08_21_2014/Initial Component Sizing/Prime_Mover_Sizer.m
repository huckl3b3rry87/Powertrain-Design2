%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%------------------------Component Sizer----------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
clear all
close all
clc


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%------------------Load the Component Data--------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
cd('Components');

% Engine
% Engine_41_kW;
% Engine_102_kW;
Engine_2rz_0410;

% Motor
Motor_75_kW;

Battery_ADVISOR;

% Vehicle
Vehicle_Parameters_4_HI_AV;
cd ..

data;                              %  Put all the data into structures
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%-------------------------Design Variables--------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
dvar.fc_trq_scale = 1;
dvar.mc_trq_scale = 1.0;
dvar.FD = 3;
dvar.G = 1;  % Also change this so that you can get to the maximum speed!

mc_max_pwr_kW =  dvar.mc_trq_scale*vinf.mc_max_pwr_kW;
dvar.module_number = ceil(4*mc_max_pwr_kW*1000*Rint_size/(Voc_size^2));
Manipulate_Data_Structure;     % Update the data based off of the new variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%--------------------------Engine Power-----------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
test = 1; % Zero Grade (Max Speed)
ttt = 1;
for FD = 1:0.05:8
    dvar.FD = FD;
    [~,~, V_max_t, ~, V_max_6,~] = Engine_Power_Sizer( param, vinf, dvar, test ); 
   
    if isempty(V_max_6);
        V_max_6 = NaN;
    end
    
    if isempty(V_max_t)
        V_max_t = NaN;
    end
    
    FD_sim(ttt) = FD;
    V_6(ttt) = V_max_6;
    V_act(ttt) =  V_max_t;
    ttt = ttt + 1;
end

[Max_Spd_temp, I_temp] = max(V_6);
FD_final = FD_sim(I_temp)*0.9;
[Junk, I] = min(abs(FD_final - FD_sim)); 
Max_Spd = V_6(I);
Final_Drive_Plot;

% Resimualte with selected gear ratio
    dvar.FD = FD_sim(I);
    [Sim_Grade,F_max_t, V_max_t, F_max_6, V_max_6, RR ] = Engine_Power_Sizer( param, vinf, dvar, test ); 
    Grade_Plot_1;

%% Now size the engine so that it meets the rest of the grade requirements
test = 2;
V_test = 55*param.mph_mps;

iii = 1;
for fc_trq_scale = 0.5:0.01:1.75
    dvar.fc_trq_scale = fc_trq_scale;
    Manipulate_Data_Structure;  % Update the data based off of the new motor power and battery module
    [Sim_Grade,F_max_t, V_max_t, F_max_6, V_max_6, RR ] = Engine_Power_Sizer( param, vinf, dvar, test );

    if isempty(V_max_t)
        V_diff(iii) = NaN;
    else
        V_diff(iii) = abs(V_max_t - V_test); % Needs to be greater!! -Change this later
    end
    fc_trq_sim(iii) =  fc_trq_scale;
    iii = iii + 1;
end
figure(3); clf
plot(fc_trq_sim,V_diff/param.mph_mps)
xlabel('fc trq scale')
xlabel('V diff')
[Min, I] = min(V_diff);
fprintf('The best fc_trq_scale is: ')
fc_trq_sim(I)  % If the fc_trq scale goes up then it will not lower the final speed

% Resimulate with selected trq scale
dvar.fc_trq_scale = fc_trq_sim(I);
Manipulate_Data_Structure;  % Update the data based off of the new motor power and battery module
[Sim_Grade,F_max_t, V_max_t, F_max_6, V_max_6, RR ] = Engine_Power_Sizer( param, vinf, dvar, test );
Grade_Plot_1;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%---------------------------Motor Gear -----------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%% Size the Motor Based off of the maximum vehicle speed - Need to find this
% based off of the new variables  and zero grade
test = 1;
n = 6;
iii = 1;

for G = 0.9:0.1:2
    clear A
    clear I
    dvar.G = G;
    [Sim_Grade,F_max_t, V_max_t, F_max_6, V_max_6, RR ] = Engine_Power_Sizer( param, vinf, dvar, test );
    V = Sim_Grade(:,n+6);
    Motor_Tractive_Effort = Sim_Grade(:,n+2);
    %     [Max_Motor_Spd ] = Max_Motor_Spd( Motor_Tractive_Effort, V );
    A =  isnan(Motor_Tractive_Effort);  % If there is a one then there is a NaN
    I = find(A==0);  % if a never = 0, all nan..
    if isempty(I)
        pause;
    end
    Max_Motor_Spd = V(I(end));
    
    V_mot = Max_Motor_Spd/param.mph_mps;
    V_eng = V_max_t/param.mph_mps;
    V_diff_mot(iii) = abs(V_mot-V_eng);
    G_sim(iii) = dvar.G;
    iii = iii + 1;
end
figure(5);clf
plot(1:(iii-1),V_diff_mot/param.mph_mps)
xlabel('itteration')
ylabel('Velocity Difference')

[Min, I] = min(V_diff_mot);
fprintf('The best G ratio is: ')
G_sim(I)

 dvar.G = G_sim(I);
 [Sim_Grade,F_max_t, V_max_t, F_max_6, V_max_6, RR ] = Engine_Power_Sizer( param, vinf, dvar, test );
 Grade_Plot;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%---------------------------Motor Size  ----------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%---Determine Acceleration Req. Based off Drive Cycles--------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
clear V_0 V_f
S = 10;                    % Number of points to calcualte acceleration at
dt_2 = 0.0002;
graph = 1;
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

dvar.mc_trq_scale = 0.4;
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





