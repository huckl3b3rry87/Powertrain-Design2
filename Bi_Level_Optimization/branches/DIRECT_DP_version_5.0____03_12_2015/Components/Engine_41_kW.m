% ADVISOR Data file:  FC_SI41_emis.M
%
% Data source: J. Dill Murrel, JDM & Associates.  
%
% Data confidence level:  
%
% Notes: 1991 Geo Metro 1.0L SI engine.
% Maximum Power 41 kW @ 5700 rpm.
% Peak Torque 81 Nm @ 3477 rpm.

% WARNING:  This data comes from transient testing on the FTP and is
% only appropriate to model transient-operation engines.
%
% Created on:  06/22/98
% By:  Tony Markel, National Renewable Energy Laboratory, Tony_Markel@nrel.gov
%
% Revision history at end of file.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPEED & TORQUE RANGES over which data is defined
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (rad/s), speed range of the engine
eng_consum_spd=[104.5 149.2 220.9 292.5 364.1 435.7 507.4 552.2 596.9]; 
eng_consum_spd_old = eng_consum_spd;
% (N*m), torque range of the engine
eng_consum_trq=[0 0.01 6.8 13.6 20.4 27.2 33.8 40.6 47.4 54.2 61 67.8 74.6 81.4]; % added zero


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUEL USE AND EMISSIONS MAPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (g/s), fuel use map indexed vertically by fc_map_spd and 
% horizontally by fc_map_trq
eng_bsfc  =[
635.7	635.7	541.4	447.2	352.9	332.2	311.4	322.4	333.5	333.5	333.5	333.5
678.4	500.1	443.8	387.4	331.1	301.8	297	283.4	269.8	358	358	358
463.4	463.4	407.6	350.1	294.3	280.8	267.3	253.9	269.8	303.2	336.7	336.7
699.1	567.9	500.3	432.7	301.4	283.9	266.3	248.7	258.8	268.8	271.9	317.9
592.9	592.9	494.6	393.4	295.1	279.4	263.6	247.9	255.2	262.5	295	322.6
667.9	524.8	381.6	351.9	322.2	304.9	287.5	270.8	290.8	310.9	330.9	330.9
630.6	630.6	522.5	411.1	303	304.4	305.8	304.2	314.5	324.8	327.7	327.7
698.4	500.5	428.6	392.7	356.8	337.9	328.4	319	328.8	338.6	333.7	333.7
751.1	637.8	521.1	407.8	393.1	378.4	363.3	348.2	318.8	340.2	340.2	340.2];
eng_bsfc =  [zeros(size(eng_bsfc,1),1), eng_bsfc(:,1), eng_bsfc];
% fuel map in g/kWh 

% (g/s), engine out HC emissions indexed vertically by fc_map_spd and
% horizontally by fc_map_trq
fc_hc_map_gpkWh =[
11.5	11.5	9.8	8.2	6.5	5.8	5.1	5.9	6.8	6.8	6.8	6.8
7.8	7	6.2	5.5	4.7	4.3	4.7	4.6	4.5	4.6	4.6	4.6
5.8	5.8	5.2	4.6	4	4	4	4	4.5	4.6	4.6	4.6
7.1	6	5.4	4.9	3.8	3.7	3.6	3.4	3.2	3	3.4	3.9
5.8	5.8	5	4.3	3.6	3.6	3.6	3.7	3.6	3.6	4	3.9
5.6	4.7	3.7	3.7	3.7	3.4	3.1	3	3.2	3.4	3.5	3.5
8.2	8.2	6.8	5.4	4.1	3.7	3.3	3.1	3.2	3.2	3.3	3.3
5.8	5.2	4.8	4.5	4.3	3.7	3.4	3.2	3.3	3.3	3.3	3.3
5.6	5.8	5.9	6.1	5.7	5.4	5	4.6	3.9	3.9	3.9	3.9];
fc_hc_map_gpkWh=[zeros(size(fc_hc_map_gpkWh,1),1), fc_hc_map_gpkWh(:,1), fc_hc_map_gpkWh];
% engine out HC in g/kWh

% (g/s), engine out CO emissions indexed vertically by fc_map_spd and
% horizontally by fc_map_trq
fc_co_map_gpkWh =[
71.8	71.8	58.8	45.7	32.7	27.5	22.4	82.9	143.3	143.3	143.3	143.3
104.4	68.3	56.8	45.3	33.9	23.3	25.7	24.2	22.8	268.6	268.6	268.6
48.3	48.3	42.9	37.2	31.8	28.6	25.4	22.2	22.8	152.9	283	283
103.1	82.7	72.2	61.8	41.5	36.9	32.3	27.8	31.1	34.4	178.5	279.9
88.1	88.1	74.2	59.9	46	41.9	37.8	33.7	34.8	35.8	158.8	264.6
96.1	74.5	52.9	51.2	49.5	45.9	42.4	34	117.9	201.8	285.7	285.7
114.6	114.6	92.1	69	46.5	60.7	74.8	129.6	195.5	261.4	277.8	277.8
60.1	63.8	64	64.1	64.3	108.2	130.2	152.2	216.7	281.2	278.1	278.1
51.8	75.2	99.3	122.8	134.9	147.1	159.7	172.2	196.6	286.6	286.6	286.6];
fc_co_map_gpkWh=[zeros(size(fc_co_map_gpkWh,1),1), fc_co_map_gpkWh(:,1), fc_co_map_gpkWh];
% engine out CO in g/kWh

% (g/s), engine out NOx emissions indexed vertically by fc_map_spd and
% horizontally by fc_map_trq
fc_nox_map_gpkWh =[
5.8	5.8	9.3	12.8	16.3	16.1	15.9	13	10.2	10.2	10.2	10.2
5.2	8.8	9.2	9.7	10.2	8	8.9	13.2	17.5	4.6	4.6	4.6
8.1	8.1	8.8	9.6	10.4	10.8	11.3	11.7	17.5	11.6	5.7	5.7
4.2	5.6	6.3	7	8.4	8.9	9.5	10.1	13.9	17.7	8.1	3.1
5.8	5.8	7.2	8.7	10.1	11	11.8	12.6	15.9	19.2	9.3	6.8
14.9	16.4	17.8	19.4	21	20.6	20.3	19.1	14.6	10.2	5.7	5.7
28.7	28.7	26.8	25	23.1	22.4	21.7	16.5	12.1	7.8	6.5	6.5
31.1	27.9	26.7	26.2	25.6	20.9	18.6	16.3	12.1	7.8	6.8	6.8
35	31.1	27.1	23.2	21.1	19.1	17	14.9	10.9	7.4	7.4	7.4];
fc_nox_map_gpkWh=[zeros(size(fc_nox_map_gpkWh,1),1), fc_nox_map_gpkWh(:,1), fc_nox_map_gpkWh];
% engine out NOx in g/kWh

% (g/s), engine out PM emissions indexed vertically by fc_map_spd and
% horizontally by fc_map_trq
fc_pm_map_gpkWh=zeros(size(eng_bsfc));

% (g/s), engine out O2 indexed vertically by fc_map_spd and
% horizontally by fc_map_trq
fc_o2_map=zeros(size(eng_bsfc));


% convert g/kWh maps to g/s maps
[T,w]=meshgrid(eng_consum_trq, eng_consum_spd);
fc_map_kW=T.*w/1000;
eng_consum_fuel=eng_bsfc.*fc_map_kW/3600;
fc_hc_map=fc_hc_map_gpkWh.*fc_map_kW/3600;
fc_co_map=fc_co_map_gpkWh.*fc_map_kW/3600;
fc_nox_map=fc_nox_map_gpkWh.*fc_map_kW/3600;
fc_pm_map=fc_pm_map_gpkWh.*fc_map_kW/3600;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LIMITS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (N*m), max torque curve of the engine indexed by fc_map_spd
eng_max_trq=[61	67.6 73.7 78.5 80.9 77.3 76.2 73.3 68.7]; 

eng_map_spd = eng_consum_spd;
[opt_vals, spd] = min(eng_bsfc,[],1);
[opt_val, optimal_trq_index] = min(opt_vals);
optimal_spd_index = spd(optimal_trq_index);
optimal_eng_trq = eng_consum_trq(optimal_trq_index);  % DP should tell us
% this
optimal_eng_spd = 3000*pi/30;

% (N*m), closed throttle torque of the engine (max torque that can be absorbed)
% indexed by fc_map_spd -- correlation from JDMA
% fc_ct_trq=4.448/3.281*(-fc_disp)*61.02/24 * ...
%    (9*(eng_consum_spd/max(eng_consum_spd)).^2 + 14 * (eng_consum_spd/max(eng_consum_spd)));
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT SCALING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% user definable mass scaling function
fc_mass_scale_fun=inline('(x(1)*fc_trq_scale+x(2))*(x(3)*fc_spd_scale+x(4))*(fc_base_mass+fc_acc_mass)+fc_fuel_mass','x','fc_spd_scale','fc_trq_scale','fc_base_mass','fc_acc_mass','fc_fuel_mass');
fc_mass_scale_coef=[1 0 1 0]; % coefficients of mass scaling function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STUFF THAT SCALES WITH TRQ & SPD SCALES (MASS AND INERTIA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fc_inertia=0.1*fc_pwr_scale;           % (kg*m^2), rotational inertia of the engine (unknown)
fc_max_pwr=(max(eng_consum_spd.*eng_max_trq)/1000); % kW     peak engine power

fc_base_mass=1.8*fc_max_pwr;            % (kg), mass of the engine block and head (base engine)
                                        %       mass penalty of 1.8 kg/kW from 1994 OTA report, Table 3 
fc_acc_mass=0.8*fc_max_pwr;             % kg    engine accy's, electrics, cntrl's - assumes mass penalty of 0.8 kg/kW (from OTA report)
fc_fuel_mass=0.6*fc_max_pwr;            % kg    mass of fuel and fuel tank (from OTA report)
fc_mass=fc_base_mass+fc_acc_mass+fc_fuel_mass; % kg  total engine/fuel system mass

fc_max_pwr_initial_kW = 41;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Define all Variables for DP
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% eng_control_trq = eng_consum_trq; % Try to define it like 0:10:max_trq later!! - See the difference
eng_control_trq  = [0:5:max(eng_consum_trq),max(eng_consum_trq)];
W_eng_min = min(eng_consum_spd);
W_eng_max = max(eng_consum_spd); 
Te_min =  min(eng_consum_trq);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Define all Variables for DP
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%


