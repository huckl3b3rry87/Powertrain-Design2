% ADVISOR Data file:  FC_INSIGHT.M
%
% Data source: 
% Fuel map from test data provided by Argonne National Laboratory (ANL) 
% Emissions maps are not included in this data file.
% 
% Data confidence level:  
%
% Notes: 
%
% Created on:  10/26/00
% By:  Ken Kelly, National Renewable Energy Laboratory, ken_kelly@nrel.gov
%
% Revision history at end of file.
% %
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % FILE ID INFO
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fc_description='Honda Insight 1.0L VTEC-E SI Engine from ANL Test Data'; 
% fc_version=2003; % version of ADVISOR for which the file was generated
% fc_proprietary=0; % 0=> non-proprietary, 1=> proprietary, do not distribute
% fc_validation=1; % 0=> no validation, 1=> data agrees with source data, 
% % 2=> data matches source data and data collection methods have been verified
% fc_fuel_type='Gasoline';
% fc_disp=1.0;  % (L), engine displacement
% fc_emis=1;      % boolean 0=no emis data; 1=emis data
% fc_cold=0;      % boolean 0=no cold data; 1=cold data exists
% disp(['Data loaded: FC_INSIGHT - ',fc_description]);
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPEED & TORQUE RANGES over which data is defined
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (rad/s), speed range of the engine 
eng_consum_spd=[800 1273 1745 2218 2691 3164 3636 4109 4582 5055 5527 6000]*2*pi/60; 

% (N*m), torque range of the engine (original data given in lb-ft and converted to N-m)
eng_consum_trq=[5.6 11.2 16.8 22.3 27.9 33.5 39.1 44.7 50.3 55.8 61.4 67.0]*1.355818; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUEL USE AND EMISSIONS MAPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fuel use map from Argonne National Laboratory
% (g/s), fuel use map indexed vertically by fc_map_spd and 
% horizontally by fc_map_trq 
eng_consum_fuel =[
  0.0962  0.1269  0.1576  0.1883  0.2191  0.2498  0.2805  0.3112  0.3610  0.4566  0.4641  0.4641 
 0.1420  0.1909  0.2398  0.2887  0.3375  0.3864  0.4353  0.4842  0.5584  0.7129  0.7383  0.7383 
 0.1871  0.2541  0.3212  0.3882  0.4552  0.5223  0.5893  0.6563  0.7533  0.9683  1.0215  1.0215 
 0.2371  0.3223  0.4075  0.4927  0.5779  0.6630  0.7482  0.8334  0.9524  1.2297  1.3207  1.3207 
 0.2953  0.3987  0.5020  0.6053  0.7087  0.8120  0.9154  1.0187  1.1591  1.5012  1.6278  1.6399 
 0.3656  0.4871  0.6086  0.7301  0.8516  0.9731  1.0946  1.2160  1.3777  1.7875  1.9363  1.9839 
 0.4521  0.5918  0.7314  0.8711  1.0107  1.1504  1.2900  1.4297  1.6124  2.0936  2.2647  2.3577 
 0.5591  0.7169  0.8747  1.0325  1.1903  1.3481  1.5059  1.6637  2.0304  2.4249  2.6182  2.7666 
 0.7038  0.8993  1.1014  1.3102  1.5255  1.7475  1.9760  2.2112  2.4530  2.7014  2.9563  3.2156 
 0.8680  1.0863  1.3123  1.5458  1.7869  2.0356  2.2919  2.5558  2.8272  3.1063  3.3930  3.5433 
 1.0663  1.3087  1.5596  1.8193  2.0877  2.3647  2.6504  2.9448  3.2479  3.5596  3.8801  3.8830 
 1.3032  1.5709  1.8484  2.1357  2.4330  2.7402  3.0572  3.3841  3.7208  4.0675  4.2459  4.2459];


% (g/s), engine out HC emissions indexed vertically by fc_map_spd and
% horizontally by fc_map_trq
fc_hc_map=zeros(size(eng_consum_fuel));
% engine out HC in g/kWh (taken from Geo 1.0L engine)

% (g/s), engine out CO emissions indexed vertically by fc_map_spd and
% horizontally by fc_map_trq 

fc_co_map=zeros(size(eng_consum_fuel));
% engine out CO in g/kWh (taken from Geo 1.0L engine)

% (g/s), engine out NOx emissions indexed vertically by fc_map_spd and
% horizontally by fc_map_trq

fc_nox_map=zeros(size(eng_consum_fuel));
% engine out NOx in g/kWh (taken from Geo 1.0L engine)

% (g/s), engine out PM emissions indexed vertically by fc_map_spd and
% horizontally by fc_map_trq
fc_pm_map=zeros(size(eng_consum_fuel));

% (g/s), engine out O2 indexed vertically by fc_map_spd and
% horizontally by fc_map_trq
fc_o2_map=zeros(size(eng_consum_fuel));

% convert g/kWh maps to g/s maps
[T,w]=meshgrid(eng_consum_trq, eng_consum_spd);
fc_map_kW=T.*w/1000;
eng_bsfc=eng_consum_fuel./fc_map_kW*3600;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LIMITS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (N*m), max torque curve of the engine indexed by fc_map_spd
% Honda Published max horsepower & torque (Honda Insight Facts Book)
%		67 hp @ 5700 rpm (50 kW @ 596.9 rad/s)
%		66 lb-ft @ 5600 4800 rpm (89.5 N-m @ 502.7 rad/s)
%
% other max torque points taken from Honda UC Davis presentation (speed vs. torque graph)
%											 
%(original data given in lb-ft and converted to N-m)
eng_max_trq=[56.9 58.2 59.5 60.7 62.0 63.2 64.5 65.7 67.0 64.3 61.5 58.6]*1.355818; 
eng_map_spd = eng_consum_spd;

% (N*m), closed throttle torque of the engine (max torque that can be absorbed)
% % indexed by fc_map_spd -- correlation from JDMA
% fc_ct_trq=4.448/3.281*(-fc_disp)*61.02/24 * ...
%    (9*(fc_map_spd/max(fc_map_spd)).^2 + 14 * (fc_map_spd/max(fc_map_spd)));

%%_scale_coef=[1 0 1 0]; % coefficients of mass scaling function
fc_pwr_scale = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STUFF THAT SCALES WITH TRQ & SPD SCALES (MASS AND INERTIA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fc_inertia=0.1*fc_pwr_scale;           % (kg*m^2), rotational inertia of the engine (unknown)
fc_max_pwr=(max(eng_consum_spd.*eng_max_trq)/1000)*fc_pwr_scale; % kW     peak engine power

fc_base_mass=60;            % (kg), SAE Automotive Engineering, Oct 99
fc_acc_mass=0.8*fc_max_pwr;             % kg    engine accy's, electrics, cntrl's - assumes mass penalty of 0.8 kg/kW (from OTA report)
fc_fuel_mass=0.6*fc_max_pwr;            % kg    mass of fuel and fuel tank (from OTA report)
fc_mass=fc_base_mass+fc_acc_mass+fc_fuel_mass; % kg  total engine/fuel system mass
fc_ext_sarea=0.3*(fc_max_pwr/100)^0.67;       % m^2    exterior surface area of engine

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Define all Variables for DP
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

optimal_eng_spd = 2000*pi/30;

% eng_control_trq = eng_consum_trq; % Try to define it like 0:10:max_trq later!! - See the difference
eng_control_trq  = [0:5:max(eng_consum_trq),eng_consum_trq(length(eng_consum_trq))];
W_eng_min = min(eng_consum_spd);
W_eng_max = max(eng_consum_spd); 
Te_min =  repmat(min(eng_consum_trq),[1 length(eng_control_trq)])';
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Define all Variables for DP
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

