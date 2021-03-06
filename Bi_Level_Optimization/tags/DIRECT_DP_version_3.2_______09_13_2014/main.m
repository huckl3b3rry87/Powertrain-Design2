clear all
close all
clc
tic 


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%----------------------------Load All Data--------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

cd('Components');
%                              ~~ Engine ~~

Engine_30_kW;
% Engine_41_kW;
% Engine_50_kW;
% Engine_73_kW;   
% Engine_95_kW;
% Engine_102_kW;
% Engine_186_kW;
% Engine_224_kW;


%                              ~~ Motor ~~
% Motor_int;
% Motor_75_kW;
% Motor_30_kW;
Motor_49_kW;
% Motor_10_kW;
% Motor_8_kW;  
% Motor_16_kW;

%                             ~~ Battery ~~
% Battery_int;  % No variation with the number of modules in this battery!!
Battery_ADVISOR;

%                              ~~ Vehicle ~~

% Vehicle_Parameters_4_HI_AV;
Vehicle_Parameters_4_HI;
% Vehicle_Parameters_8_HI_AV;
% Vehicle_Parameters_8_HI;

% Low Speed
% Vehicle_Parameters_4_low_AV;
% Vehicle_Parameters_4_low;
% Vehicle_Parameters_1_low_AV;
% Vehicle_Parameters_1_low;

cd .. 

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%-------------------Put all the data into structures----------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
data;                 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%---------------------Update the Design Variables-------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
e = 1;
dvar.FD = 4.25*e;
dvar.G = 1.7*e;
dvar.fc_trq_scale = 3.16667*e;
dvar.mc_trq_scale = 1*e;
mc_max_pwr_kW =  dvar.mc_trq_scale*vinf.mc_max_pwr_kW;
dvar.module_number = ceil(4*mc_max_pwr_kW*1000*Rint_size/(Voc_size^2));
% dvar.module_number = 12*e; 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%---------------------Update the Data-------------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
Manipulate_Data_Structure;   

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%---------------------Select Drive Cycle-------------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                              ~~ Standard ~~

cyc_name = 'HWFET';
% cyc_name = 'UDDS';
% cyc_name = 'US06';
% cyc_name = 'SHORT_CYC_HWFET';
% cyc_name = 'RAMP';
% cyc_name = 'LA92';
% cyc_name = 'CONST_65';
% cyc_name = 'CONST_45';
% cyc_name = 'COMMUTER';

% City
% cyc_name = 'INDIA_URBAN';
% cyc_name = 'MANHATTAN';
% cyc_name = 'Nuremberg';
% cyc_name = 'NYCC';
% cyc_name = 'AA_final';
%                              ~~ AV~~

% cyc_name = 'US06_AV';
% cyc_name = 'HWFET_AV';
% cyc_name = 'AA_final_AV';

[cyc_data] = Drive_Cycle(param, vinf, cyc_name );


% Check Stuff
% Engine_Plot_Test;
% eng_spd_best = optimal_eng_spd/param.rpm2rads
% break

%
RUN_TYPE = 0;  % RUN_TYPE = 1 - for DIRECT     &    RUN_TYPE = 0 - for DP only
%%
[FAIL, MPG, delta_SOC, sim] = Dynamic_Programming_func(param, vinf, dvar, cyc_data, RUN_TYPE);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%-----------------------Final Plots ect.----------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%%
if RUN_TYPE == 0
    cd('Plots')
    Main_Plot;
    Engine_Plot;
    Motor_Plot;
    %     Cost_Plot;
    Battery_Plot;
    cd ..
    MPG
    FAIL
end