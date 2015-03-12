% ADVISOR Data file:  FC_SI95.m
%
% Data source:  Reilly, Donald J.; et al. (1991). "Saturn DOHC and SOHC 
% Four Cylinder Engines." SAE paper #910676.
%
% Data confidence level:  Data in this file has been verified to represent the
% source data with acceptable accuracy.
%
% Notes:  Fuel use data was presented in Figure 4 as brake specific fuel rate
% (ug/s/W) as a function of brake mean effective pressure (x100 kPa) 
% and engine speed (x1000 rpm). From the text, the engine will provide 
% 165 Nm @ 4800 rpm.  This corresponds to 10.6x100 kPa.  Therefore, the conversion
% factor from brake mean effective pressure to torque in Nm was
% (165 Nm)/(10.6 bar)=15.6 Nm/bar.
% Maximum Power 95 kW @ 6000 rpm.
% Peak Torque 165 Nm @ 4800 rpm.

% Created on:  05/20/98
% By:  Tony Markel, National Renewable Energy Laboratory, Tony_Markel@nrel.gov
%
% Revision history at end of file.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % FILE ID INFO
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fc_description='Saturn 1.9L (95kW) DOHC SI Engine'; % one line descriptor identifying the engine
% fc_version=2003; % version of ADVISOR for which the file was generated
% fc_proprietary=0; % 0=> non-proprietary, 1=> proprietary, do not distribute
% fc_validation=0; % 1=> no validation, 1=> data agrees with source data, 
% % 2=> data matches source data and data collection methods have been verified
% fc_fuel_type='Gasoline';
% fc_disp=1.9; % (L) engine displacement
% fc_emis=0;      % boolean 0=no emis data; 1=emis data
% fc_cold=0;      % boolean 0=no cold data; 1=cold data exists
% disp(['Data loaded: FC_SI95.m - ',fc_description]);
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPEED & TORQUE RANGES over which data is defined
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (rad/s), speed range of the engine
%fc_map_spd=[eps:.5:6]*1000*2*pi/60; 
eng_consum_spd=[0.5:0.5:6]*1000*2*pi/60; 

% (N*m), torque range of the engine
%fc_map_trq=[eps:11]*15.6; 
eng_consum_trq=[1:11]*15.6; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUEL USE AND EMISSIONS MAPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (g/s), fuel use map indexed vertically by fc_map_spd and 
% horizontally by fc_map_trq
eng_bsfc=[
170   170   170   165   155   145   130   145   150   190   210   210
123   112   113   110   106   101.5 94.5  99    109   128   142   147
103   93.5  96    93    89    86.5  80    88.2  97.6  105   117   120
92    85.5  85    83.5  79    77.5  78    84.5  88    92.5  101.5 104
85    79    78    76.5  74    73    77    78.5  81.5  85    91    93
87    79.5  73    72    68.5  71    74.3  76.5  77.5  79    82    87
84    78.5  70.7  68    68    69    72    73    75.5  78.7  83    83  
90    79    68.5  68    68    68    69    70.5  75    77.5  82    82
80    80    70    69    68    68    69    73    75.8  78.6  81    81
80    80    80    73    71    72    73    74    76.5  79    80.5  80
80    80    80    80    80    80    80    80    80    80    80    80];
% fuel tbl in BSFC (ug/Ws)

% pad with speed=0 column and BMEP=0 x 100lPA column
%fc_fuel_map_ugpWs=[fc_fuel_map_ugpWs(:,1) fc_fuel_map_ugpWs];  
%fc_fuel_map_ugpWs=[fc_fuel_map_ugpWs(1,:); fc_fuel_map_ugpWs];  

% flip columns and rows the align with torque and speed vectors
eng_bsfc=eng_bsfc'; 

% convert ug/Ws to g/s and to g/kwh
[T,w]=meshgrid(eng_consum_trq, eng_consum_spd);
fc_map_kW=T.*w/1000;
eng_consum_fuel=eng_bsfc/1000000.*fc_map_kW*1000;
fc_fuel_map_gpkWh=eng_consum_fuel./fc_map_kW*3600;

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
%fc_max_trq=[7.2 7.2 8.6 9.4 10.15 10.6 10.4 10.2 10.25 10.6 10.55 10.3 9.7]*15.6; 
eng_max_trq=[7.2 8.6 9.4 10.15 10.6 10.4 10.2 10.25 10.6 10.55 10.3 9.7]*15.6; 
eng_map_spd = eng_consum_spd;

fc_pwr_scale = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STUFF THAT SCALES WITH TRQ & SPD SCALES (MASS AND INERTIA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fc_inertia=0.1*fc_pwr_scale;           % (kg*m^2), rotational inertia of the engine (unknown)
fc_max_pwr=(max(eng_consum_spd.*eng_max_trq)/1000)*fc_pwr_scale; % kW     peak engine power

fc_base_mass=1.8*fc_max_pwr;            % (kg), mass of the engine block and head (base engine)
                                        %       mass penalty of 1.8 kg/kW from 1994 OTA report, Table 3 
fc_acc_mass=0.8*fc_max_pwr;             % kg    engine accy's, electrics, cntrl's - assumes mass penalty of 0.8 kg/kW (from OTA report)
fc_fuel_mass=0.6*fc_max_pwr;            % kg    mass of fuel and fuel tank (from OTA report)
fc_mass=fc_base_mass+fc_acc_mass+fc_fuel_mass; % kg  total engine/fuel system mass
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Define all Variables for DP
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

optimal_eng_spd = 2500*pi/30;

% eng_control_trq = eng_consum_trq; % Try to define it like 0:10:max_trq later!! - See the difference
eng_control_trq  = [0:5:max(eng_consum_trq),eng_consum_trq(length(eng_consum_trq))];
W_eng_min = min(eng_consum_spd);
W_eng_max = max(eng_consum_spd); 
Te_min =  repmat(min(eng_consum_trq),[1 length(eng_control_trq)])';
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Define all Variables for DP
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
