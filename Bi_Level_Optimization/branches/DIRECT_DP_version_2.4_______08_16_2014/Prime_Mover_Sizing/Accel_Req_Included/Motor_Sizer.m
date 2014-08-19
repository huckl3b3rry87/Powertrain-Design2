clear all
close all
global m rwh ni g grade Cd Af rho Frr eng_max_trq eng_consum_spd W_eng_min W_eng_max FD G
global m_map_spd m_max_trq Wm_max Wm_min nt Voc_size Rint_size

% This is the main program that sizes everything

% First Determine The Acceleration Requirements

S = 10;  % Number of points to calcualte acceleration at
dt_2 = 0.1;

[ V_0, V_f, Acc_Final ] = Accel_Req( S, dt_2 ); % V_Final is in MPH

FD = 4.48; % Final Drive Ratio
G = 1.4;
fc_trq_scale = 1.1629;
mc_trq_scale = 1.0495;
module_number = 50;  % Do something like Adam did- updateing the mass to converge on a solution

fc_trq_scale_i = 100;  % Get things started
tt = 1; 

while  (abs(fc_trq_scale - fc_trq_scale_i)/fc_trq_scale > 0.2)
 
    fc_trq_scale_i = 1.1629;                % Set to the actual value
    
    % Next Determine the Engine Power
    
    % Select a vehicle and components
    cd('Components');
    Engine_41_kW;
    % Engine_102_kW;
    Battery_ADVISOR;
    Motor_75_kW;
    Vehicle_Parameters_4_HI_AV;
    % Vehicle_Parameters_4_HI;
    cd ..
    
    % Desired Grade and Top Speed Performance - (Engine) - For Low Grades
    % Test 1
    V_eng = 50*mph_mps;
    alpha_eng = 5*pi/180;
    Pbase1 = Engine_Base_Power(V_eng,alpha_eng);
    
    % Test 2
    V_eng = 75*mph_mps;
    alpha_eng = 0*pi/180;
    Pbase2 = Engine_Base_Power(V_eng, alpha_eng);
    
    % Initial Enigne Power
    Pbase_kW =  max([Pbase1;Pbase2]);
    
    % Calcualte the fc_trq_scale
    Paux_kW = Paux/1000;
    fc_trq_scale = Pbase_kW/fc_max_pwr + Paux_kW/fc_max_pwr_initial_kW;  % For now neglect the effect that changing the mass of the engine has on this
    
    % Size the Motor based off the acceleration requirements
    
    % V_0 = 1;
    % V_f = 2;
    % S = 2;
    x_sim_full =[];
    V_sim_full = [];
    acc_sim_full = [];
    Teng_sim_full =[];
    Tm_sim_full =[];
    We_sim_full= [];
    Wm_sim_full= [];
    time_sim_full = [];
    index(1) = 1;
    
    for i = 1:(S-1)    % Pass the Velocities in MPH
        i
        [ mc_trq_scale, Module_Number,Sim_Variables] = Motor_Scaling_Function(V_0(i), V_f(i),fc_trq_scale, dt_2);
        mc_trq_scale_save(i) = mc_trq_scale;
        Module_Number_Save(i) = Module_Number;
        x_sim_full = [x_sim_full; Sim_Variables(2:end,1)];
        V_sim_full = [V_sim_full; Sim_Variables(2:end,2)];
        acc_sim_full = [acc_sim_full; Sim_Variables(2:end,3)];
        Teng_sim_full = [Teng_sim_full; Sim_Variables(2:end,4)];
        Tm_sim_full =  [Tm_sim_full; Sim_Variables(2:end,5)];
        We_sim_full=  [We_sim_full; Sim_Variables(2:end,6)];
        Wm_sim_full= [Wm_sim_full; Sim_Variables(2:end,7)];
        time_sim_full = [time_sim_full; Sim_Variables(2:end,8)];
        
        index(i+1) = [index(i) + length(Sim_Variables(2:end,1))];
    end
    %%
    e = 6;
    start = index(e);
    stop = index(e+1)-1;
    
    x_sim = x_sim_full(start:stop);
    V_sim =  V_sim_full(start:stop);
    acc_sim = acc_sim_full(start:stop);
    Teng_sim =  Teng_sim_full(start:stop);
    Tm_sim = Tm_sim_full(start:stop);
    We_sim = We_sim_full(start:stop);
    Wm_sim =  Wm_sim_full(start:stop);
    time_sim = time_sim_full(start:stop);
    
    Kinematics_Plot;
    Kinetics_Plot;
    
    mc_trq_scale =  max(mc_trq_scale_save)
    module_number =  max(Module_Number_Save)
    abs(fc_trq_scale - fc_trq_scale_i)/fc_trq_scale
    fc
    figure(3);
    plot(1:S-1,mc_trq_scale_save)
    ylabel('mc trq scale')
    xlabel('Test Runs')
    tt = tt + 1
    pause;
end












