% ADVISOR Data file:  FC_SI186.M
%
% Data source:  Fuel consumption data was provided by Keith Vertin (NREL).
% Performance data from John Deere spec sheet.
%
% Data confidence level:  no comparison has been performed
%
% Notes: Turbocharged spark ignited compressed natural gas engine.
% Max torque - 1085 Nm @ 1300 rpm
% Max power - 186 kW @ 2200 rpm
% Idle speed - 650 rpm
% Bore X Stroke - 116 mm X 129 mm
%
% Created on:  05/20/98
% By:  Tony Markel, National Renewable Energy Laboratory, Tony_Markel@nrel.gov
%
% Revision history at end of file.

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % FILE ID INFO
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fc_description='John Deere 8.1L (186kW) CNG SI Engine'; % one line descriptor identifying the engine
% fc_version=2003; % version of ADVISOR for which the file was generated
% fc_proprietary=0; % 0=> non-proprietary, 1=> proprietary, do not distribute
% fc_validation=0; % 0=> no validation, 1=> data agrees with source data, 
% % 2=> data matches source data and data collection methods have been verified
% fc_fuel_type='CNG';
% fc_disp=8.1; % (L) engine displacement
% fc_emis=0;      % boolean 0=no emis data; 1=emis data
% fc_cold=0;      % boolean 0=no cold data; 1=cold data exists
% disp(['Data loaded: FC_SI186.M - ',fc_description]);
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPEED & TORQUE RANGES over which data is defined
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (rad/s), speed range of the engine
eng_consum_spd=[600 800 1000 1200 1400 1600 1800 2000 2200]*2*pi/60; 

% (N*m), torque range of the engine
eng_consum_trq=[1 50 100 200 300 400 500 600 700 800]*1055/778.16; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUEL USE AND EMISSIONS MAPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (g/s), fuel use map indexed vertically by fc_map_spd and 
% horizontally by fc_map_trq
fc_map_bte = [
8	13 21	28.5	29.5	35.5	35.5	36	36.5	36.5
10	15 20	30	35	37.3	37.6	37.8	37.6	37.3
13	17 22	31	36	38	38.2	38.5	38.5	38.5
14	20 26.5	30.5	35.6	37.1	37.5	37.8	38.2	39.1
16	20 23.5	33.2	35.5	37	37.8	38.5	39.2	39.5
18	20.5 23	32.4	35.1	37.2	37.9	38.5	39.2	39.5
18.5 21	23.9	31.7	34.9	36.6	37.6	38.2	38.9	39.3
19	21.2 24	31	34.5	35.8	37.3	37.9	38.3 38.3
19	20.5 21	30	33.6	36.1	37	36.9	37.5	37.5]./100; 
% brake thermal efficiency in decimal percent

% convert map from bte to gps
[T,w]=meshgrid(eng_consum_trq,eng_consum_spd);
fc_map_kW=T.*w/1000;
fc_fuel_lhv=47.0*1000; % (J/g), lower heating value of CNG
eng_consum_fuel=1/fc_fuel_lhv./fc_map_bte.*fc_map_kW*1000;
eng_bsfc=eng_consum_fuel./fc_map_kW*3600;


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
eng_max_trq=[400 500 687.5 785 790 762.5 715 665 597]*1055/778.16;
% unknown
fc_pwr_scale = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STUFF THAT SCALES WITH TRQ & SPD SCALES (MASS AND INERTIA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fc_inertia=0.1*fc_pwr_scale;   % (kg*m^2),  rotational inertia of the engine (unknown)
fc_max_pwr=(max(eng_consum_spd.*eng_max_trq)/1000)*fc_pwr_scale; % kW     peak engine power

fc_base_mass=2.8*fc_max_pwr;            % (kg), mass of the engine block and head (base engine)
                                        %       assuming a mass penalty of 1.8 kg/kW from S. Sluder (ORNL) estimate of 300 lb 
fc_acc_mass=0.8*fc_max_pwr;             % kg    engine accy's, electrics, cntrl's - assumes mass penalty of 0.8 kg/kW (from 1994 OTA report, Table 3)
fc_fuel_mass=0.6*fc_max_pwr;            % kg    mass of fuel and fuel tank (from 1994 OTA report, Table 3)
fc_mass=fc_base_mass+fc_acc_mass+fc_fuel_mass; % kg  total engine/fuel system mass
%fc_mass=667;                            % (kg), mass of the engine
fc_ext_sarea=0.3*(fc_max_pwr/100)^0.67;       % m^2    exterior surface area of engine

fc_inertia=fc_mass*(1/3+1/3*2/3)*(0.15^2) ;
% assumes 1/3 purely rotating mass, 1/3 purely oscillating, and 1/3 stationary
% and crank radius of 0.15m, 2/3 of oscilatting mass included in rotational inertia calc
% correlation from Bosch handbook p.379
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Define all Variables for DP
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
eng_map_spd = eng_consum_spd;
optimal_eng_spd = 1400*pi/30;

% eng_control_trq = eng_consum_trq; % Try to define it like 0:10:max_trq later!! - See the difference
eng_control_trq  = [0:5:max(eng_consum_trq),eng_consum_trq(length(eng_consum_trq))];
W_eng_min = min(eng_consum_spd);
W_eng_max = max(eng_consum_spd); 
Te_min =  repmat(min(eng_consum_trq),[1 length(eng_control_trq)])';
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Define all Variables for DP
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
