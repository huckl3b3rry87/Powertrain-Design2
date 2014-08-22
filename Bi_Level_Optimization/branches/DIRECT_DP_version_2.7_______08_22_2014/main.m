clear all
close all
clc
tic 

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
% Motor_49_kW;

%                             ~~ Battery ~~
% Battery_int;  % No variation with the number of modules in this battery!!
Battery_ADVISOR;

%                              ~~ Vehicle ~~

Vehicle_Parameters_4_HI_AV;
% Vehicle_Parameters_4_HI;
% Vehicle_Parameters_8_HI_AV;
% Vehicle_Parameters_8_HI;

cd .. 

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%-------------------Put all the data into structures----------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
data;                 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%---------------------Update the Design Variables-------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
dvar.FD = 3.65;
dvar.G = 1.6;
dvar.fc_trq_scale =  0.87;
dvar.mc_trq_scale = 1;
mc_max_pwr_kW =  dvar.mc_trq_scale*vinf.mc_max_pwr_kW;
% dvar.module_number = ceil(4*mc_max_pwr_kW*1000*Rint_size/(Voc_size^2));
dvar.module_number = 19;  
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%---------------------Update the Data-------------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
Manipulate_Data_Structure;   

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%---------------------Select Drive Cycle-------------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% 
% cyc_name = 'HWFET';
% cyc_name = 'UDDS';
% cyc_name = 'US06';
% cyc_name = 'SHORT_CYC_HWFET';
% cyc_name = 'RAMP';
% cyc_name = 'LA92';
% cyc_name = 'CONST_65';
% cyc_name = 'CONST_45';
% cyc_name = 'COMMUTER';

cyc_name = 'HWFET_AV';

[cyc_data] = Drive_Cycle(param, vinf, cyc_name );

RUN_TYPE = 0;  % RUN_TYPE = 1 - for DIRECT     &    RUN_TYPE = 0 - for DP only

[FAIL, MPG, delta_SOC, sim] = Dynamic_Programming_func(param, vinf, dvar, cyc_data, RUN_TYPE);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%-----------------------Final Plots ect.----------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

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