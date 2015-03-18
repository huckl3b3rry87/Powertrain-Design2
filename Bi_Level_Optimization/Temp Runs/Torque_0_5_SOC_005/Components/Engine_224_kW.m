% ADVISOR Data file:  FC_SI224.m
%
% Data source:  Cummins Engine Company published spec sheets
%
% Data confidence level:  no comparison has been performed
%
% Notes:  Data provided included, mass, max torque vs. spd, and fuel consumption @ full load vs. spd. 
% All other information has been estimated!!
%
% Created on:  05/28/98
% By:  Tony Markel, National Renewable Energy Laboratory, Tony_Markel@nrel.gov
%
% Revision history at end of file.

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % FILE ID INFO
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fc_description='Cummins L10-300G (224 kW) LNG Engine'; % one line descriptor identifying the engine
% fc_version=2003; % version of ADVISOR for which the file was generated
% fc_proprietary=0; % 0=> non-proprietary, 1=> proprietary, do not distribute
% fc_validation=0; % 0=> no validation, 1=> data agrees with source data, 
% % 2=> data matches source data and data collection methods have been verified
% fc_fuel_type='LNG';
% fc_disp=10.0; % (L) engine displacement
% fc_emis=0;      % boolean 0=no emis data; 1=emis data
% fc_cold=0;      % boolean 0=no cold data; 1=cold data exists
% disp(['Data loaded: FC_SI224.M - ',fc_description]);
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPEED & TORQUE RANGES over which data is defined
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (rad/s), speed range of the engine
fc_map_spd=[1100 1200 1300 1400 1500 1600 1700 1800 1900 2100]*2*pi/60; 

% (N*m), torque range of the engine
fc_map_trq=[100:100:1300]; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUEL USE AND EMISSIONS MAPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (g/s), fuel use map indexed vertically by fc_map_spd and 
% horizontally by fc_map_trq
fc_fuel_map =[
    1.2852    1.9885    2.5767    3.0856    3.5948    4.1828    4.8270    5.4884    6.1430    6.7992    7.4484    8.1093    8.7676
    1.2877    2.0533    2.8080    3.3731    3.9572    4.6283    5.3513    6.0777    6.8001    7.5045    8.1638    8.8882    9.6096
    1.3966    2.2200    2.8613    3.5247    4.1741    4.8774    5.6139    6.3528    7.0863    7.8028    8.5033    9.1846    9.9301
    1.5686    2.4833    3.0103    3.8167    4.6021    5.3766    6.1786    6.9683    7.7448    8.5027    9.2707   10.0933   10.9126
    1.7091    2.7450    3.3556    4.2119    5.0548    5.8836    6.7574    7.6280    8.4837    9.3145   10.1599   11.0614   11.9593
    1.8879    3.0738    3.7808    4.7087    5.6232    6.5198    7.4894    8.4621    9.4181   10.3392   11.2828   12.2840   13.2811
    2.0033    3.2960    4.1255    5.1018    6.0886    7.0659    8.1080    9.1487   10.1827   11.1799   12.1891   13.2706   14.3478
    2.0802    3.4543    4.3861    5.4080    6.4462    7.5027    8.5876    9.6743   10.7684   11.8236   12.9799   14.1317   15.2788
    2.1574    3.6122    4.6014    5.6756    6.7591    7.8925    8.9959   10.1042   11.2495   12.3796   13.5904   14.7963   15.9974
    2.4188    4.2045    5.2453    6.3815    7.5129    8.7041    9.8927   11.0890   12.4006   13.7025   15.0426   16.3774   17.7067];
% The above map was generate using the lines of code encased in an 
% if exist('fc_fuel_map') statement provided at the end of this file.
% The map takes on the shape of the John Deere 8.1L NG engine while maintaining 
% absolute fuel consumption data along its max trq curve.

% generate the BSFC map
[T,w]=meshgrid(fc_map_trq,fc_map_spd);
fc_map_kW=T.*w/1000;
fc_fuel_map_gpkWh=fc_fuel_map./fc_map_kW*3600;

% pad maps with idle conditions
fc_idle_spd=650;
fc_map_spd=[fc_idle_spd*pi/30 fc_map_spd];
fc_fuel_map=[fc_fuel_map(1,:)/1.15; fc_fuel_map];
fc_fuel_map_gpkWh=[fc_fuel_map_gpkWh(1,:)*1.5; fc_fuel_map_gpkWh];

% pad maps with overspd condition
%fc_max_overspd=2300;
%fc_map_spd=[fc_map_spd fc_max_overspd*pi/30];
%fc_fuel_map=[fc_fuel_map; fc_fuel_map(end,:)*1.15];
%fc_fuel_map_gpkWh=[fc_fuel_map_gpkWh; fc_fuel_map_gpkWh(end,:)*1.5];

% pad maps with a zero spd and zero torque column
%fc_map_spd=[eps fc_map_spd];
%fc_map_trq=[eps fc_map_trq];
%fc_fuel_map=[fc_fuel_map(:,1)/2 fc_fuel_map];
%fc_fuel_map=[fc_fuel_map(1,:)/2; fc_fuel_map];
%fc_fuel_map_gpkWh=[fc_fuel_map_gpkWh(:,1)*1.5 fc_fuel_map_gpkWh];
%fc_fuel_map_gpkWh=[fc_fuel_map_gpkWh(1,:)*1.5; fc_fuel_map_gpkWh];


% (g/s), engine out HC emissions indexed vertically by fc_map_spd and
% horizontally by fc_map_trq
fc_hc_map=zeros(size(fc_fuel_map));

% (g/s), engine out HC emissions indexed vertically by fc_map_spd and
% horizontally by fc_map_trq
fc_co_map=zeros(size(fc_fuel_map));

% (g/s), engine out HC emissions indexed vertically by fc_map_spd and
% horizontally by fc_map_trq
fc_nox_map=zeros(size(fc_fuel_map)); 

% (g/s), engine out PM emissions indexed vertically by fc_map_spd and
% horizontally by fc_map_trq
fc_pm_map=zeros(size(fc_fuel_map));

% (g/s), engine out O2 indexed vertically by fc_map_spd and
% horizontally by fc_map_trq
fc_o2_map=zeros(size(fc_fuel_map));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LIMITS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (N*m), max torque curve of the engine indexed by fc_map_spd
fc_max_trq=[1139 1180 1220 1191 1163 1134 1106 1077 1047 1017]; 
% pad for idle condition and overspd condition
%fc_max_trq=[fc_max_trq(1)*0.95 fc_max_trq fc_max_trq(end)*0.9]; 
% pad for idle condition
fc_max_trq=[fc_max_trq(1)*0.95 fc_max_trq]; 

fc_pwr_scale = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STUFF THAT SCALES WITH TRQ & SPD SCALES (MASS AND INERTIA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fc_inertia=0.1*fc_pwr_scale;   % (kg*m^2),  rotational inertia of the engine (unknown)
% inertia must be greater than 5 with the specified accessory loads (tm:5/27/99)
%fc_inertia=10*fc_pwr_scale;   % (kg*m^2),  rotational inertia of the engine (unknown)
fc_max_pwr=(max(fc_map_spd.*fc_max_trq)/1000)*fc_pwr_scale; % kW     peak engine power

fc_base_mass=2.8*fc_max_pwr;            % (kg), mass of the engine block and head (base engine)
                                        %       assuming a mass penalty of 1.8 kg/kW from S. Sluder (ORNL) estimate of 300 lb 
fc_acc_mass=0.8*fc_max_pwr;             % kg    engine accy's, electrics, cntrl's - assumes mass penalty of 0.8 kg/kW (from 1994 OTA report, Table 3)
fc_fuel_mass=0.6*fc_max_pwr;            % kg    mass of fuel and fuel tank (from 1994 OTA report, Table 3)
fc_mass=fc_base_mass+fc_acc_mass+fc_fuel_mass; % kg  total engine/fuel system mass
%fc_mass=983; % (kg), mass of the engine (reported mass)
% fc_ext_sarea=0.3*(fc_max_pwr/100)^0.67;       % m^2    exterior surface area of engine
% 
% fc_inertia=fc_mass*(1/3+1/3*2/3)*(0.15^2); 
% assumes 1/3 purely rotating mass, 1/3 purely oscillating, and 1/3 stationary
% and crank radius of 0.25m, 2/3 of oscilatting mass included in rotational inertia calc
% correlation from Bosch handbook p.37

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Define all Variables for DP
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
[opt_vals, spd] = min(eng_bsfc,[],1);
[opt_val, optimal_trq_index] = min(opt_vals);
optimal_spd_index = spd(optimal_trq_index);
optimal_eng_trq = eng_consum_trq(optimal_trq_index);  % DP should tell us this
optimal_eng_spd = eng_consum_spd(optimal_spd_index);

% eng_control_trq = eng_consum_trq; % Try to define it like 0:10:max_trq later!! - See the difference
eng_control_trq  = [0:5:max(eng_consum_trq),eng_consum_trq(length(eng_consum_trq))];
W_eng_min = min(eng_consum_spd);
W_eng_max = max(eng_consum_spd); 
Te_min =  repmat(min(eng_consum_trq),[1 length(eng_control_trq)])';
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Define all Variables for DP
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
