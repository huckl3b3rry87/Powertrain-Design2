%%%%% ADVISOR data file:  MC_AC30
%
% Data source:
% Siemens electric motor 30 kW
% motor/controller combination at 216 V, dated 1996
%
% Data confidence level:
%
% Notes:
% Efficiency/loss data appropriate to a 216-V system.
%
% Created on:  24/02/99
% By: G. Villosio  
%
% Revision history at end of file.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TORQUE AND SPEED ranges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N-m, torque vector corresponding to columns of efficiency & loss maps

mc_map_trq=[-125 -105 -95 -85 -75 -65 -55 -45 -30 -15 -5 ...
	0 5 15 30 45 55 65 75 85 95 105 125];

% rad/s, speed vector corresponding to rows of efficiency & loss maps
mc_map_spd=[0 1500 2200 3000 4500 6000 7500 9000]*(2*pi)/60;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOSSES AND EFFICIENCIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% motor and inverter efficiency
%-125 -105  -95   -85   -75   -65   -55   -45   -30   -15   -5    0     5     15    30    45    55    65    75    85    95    105   125
mc_eff_map =[...
0.32	0.37	0.4	0.46	0.49	0.51	0.50	0.49	0.48	0.45	0.40	0.1	0.40	0.45	0.48	0.49	0.50	0.51	0.49	0.46	0.40	0.37	0.32    % 0
0.707	0.71	0.74	0.74	0.75	0.754	0.76	0.755	0.75	0.72	0.62	0.1	0.62	0.72	0.75	0.755	0.76	0.754	0.75	0.74	0.74	0.71	0.707   % 1500
0.707	0.71	0.74	0.74	0.75	0.754	0.76	0.755	0.75	0.72	0.62	0.1	0.62	0.72	0.75	0.755	0.76	0.754	0.75	0.74	0.74	0.71	0.707   % 2200
0.707	0.82	0.83	0.834	0.844	0.842	0.842	0.843	0.83	0.785	0.70	0.1	0.70	0.785	0.83	0.843	0.842	0.842	0.844	0.834	0.83	0.822	0.707   % 3000
0.707	0.82	0.83	0.834	0.844	0.803	0.84	0.847	0.855	0.817	0.71	0.1	0.71	0.817	0.855	0.847	0.84	0.803	0.844	0.834	0.83	0.822	0.707   % 4500
0.707	0.82	0.83	0.834	0.844	0.803	0.84	0.741	0.84	0.84	0.74	0.1	0.74	0.84	0.84	0.741 0.84	0.803	0.844	0.834 0.83	0.822	0.707   % 6000
0.707	0.82	0.83	0.834	0.844	0.803	0.84	0.741	0.825	0.825	0.76	0.1	0.76	0.825	0.825	0.741	0.84	0.803	0.844	0.834	0.83	0.822	0.707   % 7500
0.707	0.82	0.83	0.834	0.844	0.803	0.84	0.741	0.825	0.825	0.76	0.1	0.76	0.825	0.825	0.741	0.84	0.803	0.844	0.834	0.83	0.822	0.707]; % 9000
																				
% CONVERT EFFICIENCY MAP TO INPUT POWER MAP
%% find indices of well-defined efficiencies (where speed and torque > 0)
pos_trqs=find(mc_map_trq>0);
pos_spds=find(mc_map_spd>0);

%% compute losses in well-defined efficiency area
[T1,w1]=meshgrid(mc_map_trq(pos_trqs),mc_map_spd(pos_spds));
mc_outpwr1_map=T1.*w1;
mc_losspwr_map=(1./mc_eff_map(pos_spds,pos_trqs)-1).*mc_outpwr1_map;

%% to compute losses in entire operating range
%% ASSUME that losses are symmetric about zero-torque axis, and
%% ASSUME that losses at zero torque are the same as those at the lowest positive
%% torque, and
%% ASSUME that losses at zero speed are the same as those at the lowest positive
%% speed
mc_losspwr_map=[fliplr(mc_losspwr_map) mc_losspwr_map(:,1) mc_losspwr_map];
mc_losspwr_map=[mc_losspwr_map(1,:);mc_losspwr_map];

%% compute input power (power req'd at electrical side of motor/inverter set)
[T,w]=meshgrid(mc_map_trq,mc_map_spd);
mc_outpwr_map=T.*w;
mc_inpwr_map=mc_outpwr_map+mc_losspwr_map;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LIMITS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% max torque curve of the motor indexed by mc_map_spd
%           0   1500 2200 3000 4500  6000 7500 9000 
mc_max_trq=[124 124  124  92   62    45   35   28]; % (N*m)

mc_max_gen_trq=-1*[124 124  124  92   62    45   35   28]; % (N*m), estimate


% maximum overtorque capability (not continuous, because the motor would overheat) 
mc_overtrq_factor=1; % (--), estimated

mc_max_crrnt=280; % (A), maximum current allowed by the controller and motor
mc_min_volts=140; % (V), minimum voltage allowed by the controller and motor

mc_max_pwr_kW = max(mc_max_trq.*mc_map_spd)/1000;

mc_max_pwr_initial_kW = 30;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT SCALING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % (--), used to scale mc_map_spd to simulate a faster or slower running motor 
% mc_spd_scale=1.0;
% 
% % (--), used to scale mc_map_trq to simulate a higher or lower torque motor
% mc_trq_scale=1.0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OTHER DATA		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%													
% mc_inertia=0.24;		% rotor's rotational inertia, kg-m^2

% Weight of motor power/sum of motor and controller mass																		
mc_mass=41 + 17; % 

% motor/controller thermal model
mc_th_calc=1;                    % --     0=no mc thermal calculations, 1=do calcs
mc_cp=430;                       % J/kgK  ave heat capacity of motor/controller (estimate: ave of SS & Cu)
mc_tstat=45;                     % C      thermostat temp of motor/controler when cooling pump comes on
mc_area_scale=(mc_mass/91)^0.7;  % --     if motor dimensions are unknown, assume rectang shape and scale vs AC75
mc_sarea=0.4*mc_area_scale;      % m^2    total module surface area exposed to cooling fluid (typ rectang module)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Define all Variables for DP
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
m_map_spd = mc_map_spd;
m_map_trq = mc_map_trq;

m_max_trq = mc_max_trq;
m_max_gen_trq = mc_max_gen_trq;

m_eff_map = mc_eff_map;

Wm_min = -max(mc_map_spd);
Wm_max =  max(mc_map_spd);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Define all Variables for DP
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
