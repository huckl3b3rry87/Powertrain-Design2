% ADVISOR Data file:  FC_SI30.M
%
% Data source:  Data provided by University of Maryland for their 
% converted 3-cylinder Subaru Justy engine.
%
% Data confidence level:  
%
% Notes:  Engine fuel injectors and cylinder heads were converted by 
% University of Maryland to use M85 fuel. Data was provided as efficiency 
% as a function of torque (ft-lb) and speed (rpm). It is converted to a 
% fuel use table, g/kWh as a function of torque (Nm) and speed (rad/s).
% This is the same engine as that modeled in a_geo1l.m  the specific power here is 
% the value for that model scaled by the ratio of the respective maximum powers.
% Maximum Power 30 kW @ 3500 rpm.
% Fuel Converter Peak Torque 81 Nm @ 3500 rpm.
%
% Created on:  06/17/98
% By:  Tony Markel, National Renewable Energy Laboratory, Tony_Markel@nrel.gov
%
% Revision history at end of file.
% 
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILE ID INFO
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fc_description='Univ. of Maryland 30kW M85 Engine'; % one line descriptor identifying the engine
% fc_version=2003; % version of ADVISOR for which the file was generated
% fc_proprietary=0; % 0=> non-proprietary, 1=> proprietary, do not distribute
% fc_validation=0; % 1=> no validation, 1=> data agrees with source data, 
% 2=> data matches source data and data collection methods have been verified
% fc_fuel_type='M85';
% fc_disp=1.0; % (L) engine displacement
% fc_emis=0;      % boolean 0=no emis data; 1=emis data
% fc_cold=0;      % boolean 0=no cold data; 1=cold data exists
% disp(['Data loaded: FC_SI30.m - ',fc_description]);
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPEED & TORQUE RANGES over which data is defined
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (rad/s), speed range of the engine
% fc_map_spd=[0,1000:500:4500]*pi/30; 
eng_consum_spd=[1000:500:4500]*pi/30; 
eng_map_spd=eng_consum_spd;
% (N*m), torque range of the engine
eng_consum_trq=[5:5:60]*1.356; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUEL USE AND EMISSIONS MAPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (g/s), fuel use map indexed vertically by fc_map_spd and 
% horizontally by fc_map_trq
fc_eff_map=[0.118	0.112	0.107	0.118	0.113	0.102	0.095	0.084
0.165	0.165	0.167	0.168	0.168	0.158	0.151	0.148
0.199	0.199	0.206	0.211	0.213	0.199	0.191	0.179
0.229	0.229	0.238	0.239	0.236	0.228	0.211	0.190
0.253	0.253	0.251	0.250	0.252	0.244	0.226	0.209
0.265	0.265	0.266	0.268	0.267	0.261	0.247	0.228
0.275	0.275	0.279	0.282	0.275	0.273	0.258	0.232
0.284	0.284	0.290	0.293	0.287	0.282	0.257	0.224
0.290	0.290	0.294	0.298	0.297	0.289	0.256	0.224
0.294	0.294	0.302	0.307	0.305	0.294	0.256	0.224
0.294	0.294	0.302	0.312	0.311	0.306	0.256	0.224
0.294	0.294	0.302	0.312	0.311	0.311	0.256	0.224];

% fuel converter efficiency map in decimal percent
fc_eff_map=fc_eff_map'; % invert map 


% convert efficiency to g/s fuel consumption
% specific calorific values in MJ/kg from Bosch handbook 3rd ed. page 232
meth_spec_cal=19.7;                                % specific calorific value for methanol, MJ/kg
gas_spec_cal=42.7;                                 % specific calorific value for gasoline, MJ/kg
M85_spec_cal=.85*meth_spec_cal+.15*gas_spec_cal;   % specific calorific value for M85, MJ/kg
fc_fuel_lhv= 1000*M85_spec_cal;                    % (J/g), lower heating value of the fuel 23150 J/g
[T,w]=meshgrid(eng_consum_trq, eng_consum_spd);
fc_map_kW=(T.*w/1000);
eng_consum_fuel=(fc_map_kW./fc_eff_map)*(1000/fc_fuel_lhv);
eng_bsfc=eng_consum_fuel./fc_map_kW*3600;

% minimum speed in above data is 1000 rpm--assume BSFC at zero speed is 4*BSFC at 1000 rpm
%fc_fuel_map=[1*fc_fuel_map(1,:) ; fc_fuel_map];

% minimum torque in above data is 6.8 Nm--assume BSFC at zero torque is 4*BSFC at 6.8 Nm
%fc_fuel_map=[1*fc_fuel_map(:,1) fc_fuel_map];

% (g/s), engine out HC emissions indexed vertically by fc_map_spd and
% horizontally by fc_map_trq
fc_hc_map=zeros(size(eng_consum_fuel));

% (g/s), engine out HC emissions indexed vertically by fc_map_spd and
% horizontally by fc_map_trq
fc_co_map=zeros(size(eng_consum_fuel));

% (g/s), engine out HC emissions indexed vertically by fc_map_spd and
% horizontally by fc_map_trq
fc_nox_map=zeros(size(eng_consum_fuel)); 

% (g/s), engine out PM emissions indexed vertically by fc_map_spd and
% horizontally by fc_map_trq
fc_pm_map=zeros(size(eng_consum_fuel));

% (g/s), engine out O2 indexed vertically by fc_map_spd and
% horizontally by fc_map_trq
fc_o2_map=zeros(size(eng_consum_fuel));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LIMITS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (N*m), max torque curve of the engine indexed by fc_map_spd
%fc_max_trq=[0 5 50 50 55 55 60 45 40]*1.356; 
eng_max_trq=[5 50 50 55 55 60 45 40]*1.356; 

fc_pwr_scale = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STUFF THAT SCALES WITH TRQ & SPD SCALES (MASS AND INERTIA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fc_inertia=0.1*fc_pwr_scale;           % (kg*m^2), rotational inertia of the engine (unknown)
fc_max_pwr=(max(eng_consum_spd.*eng_max_trq)/1000)*fc_pwr_scale; % kW     peak engine power

fc_base_mass=1.38*1.8*fc_max_pwr;            % (kg), mass of the engine block and head (base engine)
                                        %       mass penalty of 1.8 kg/kW from 1994 OTA report, Table 3 
fc_acc_mass=1.38*0.8*fc_max_pwr;             % kg    engine accy's, electrics, cntrl's - assumes mass penalty of 0.8 kg/kW (from OTA report)
fc_fuel_mass=1.38*0.6*fc_max_pwr;            % kg    mass of fuel and fuel tank (from OTA report)
fc_mass=fc_base_mass+fc_acc_mass+fc_fuel_mass; % kg  total engine/fuel system mass
% estimated assuming a specific energy of 0.285 kW/kg (including exhaust and
% other accessories) and a 0.726 ratio of M85 max power to gasoline max power
%fc_mass=max(fc_map_spd.*fc_max_trq)/(0.726*0.285*1000); % (kg)
fc_ext_sarea=0.3*(fc_max_pwr/100)^0.67;       % m^2    exterior surface area of engine

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Define all Variables for DP
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

optimal_eng_spd = eng_consum_spd(5);

% eng_control_trq = eng_consum_trq; % Try to define it like 0:10:max_trq later!! - See the difference
eng_control_trq  = [0:5:max(eng_consum_trq),eng_consum_trq(length(eng_consum_trq))];
W_eng_min = min(eng_consum_spd);
W_eng_max = max(eng_consum_spd); 
Te_min =  repmat(min(eng_consum_trq),[1 length(eng_control_trq)])';
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Define all Variables for DP
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
