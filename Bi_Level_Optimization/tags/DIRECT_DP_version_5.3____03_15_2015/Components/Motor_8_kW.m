%%%%% ADVISOR data file:  MC_PM8
%
% Data source:
% Unique Mobility specification sheet dated 10/1/94
%
% Data confirmation:
%
% Caveats:
% Efficiency/loss data appropriate to a 100-V system.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TORQUE AND SPEED ranges at which data is available
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N-m, torque vector corresponding to columns of efficiency & loss maps
mc_map_trq=[-180 -160 -140 -120 -100 -80 -60 -40 -20 ...
	0 20 40 60 80 100 120 140 160 180]*4.448/3.281/12;

% rad/s, speed vector corresponding to rows of efficiency & loss maps
mc_map_spd=[0 500 1000 1500 2000 2500 3000 3500 4000 4500 5000]*(2*pi)/60;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOSSES AND EFFICIENCIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mc_eff_map =[...
0.1	0.1	0.1	0.1	0.1	0.1	0.1	0.1	0.1	0.1	0.1	0.1	0.1	0.1	0.1	0.1	0.1	0.1	0.1
0.18	0.28	0.41	0.52	0.42	0.36	0.38	0.3	0.2	0.1	0.2	0.3	0.38	0.36	0.42	0.52	0.41	0.28	0.18
0.46	0.57	0.63	0.71	0.7	0.69	0.71	0.6	0.34	0.1	0.34	0.6	0.71	0.69	0.7	0.71	0.63	0.57	0.46
0.63	0.71	0.74	0.79	0.8	0.82	0.85	0.76	0.48	0.1	0.48	0.76	0.85	0.82	0.8	0.79	0.74	0.71	0.63
0.74	0.8	0.82	0.835	0.85	0.855	0.865	0.855	0.61	0.1	0.61	0.855	0.865	0.855	0.85	0.835	0.82	0.8	0.74
0.81	0.835	0.85	0.855	0.865	0.87	0.875	0.86	0.72	0.1	0.72	0.86	0.875	0.87	0.865	0.855	0.85	0.835	0.81
0.84	0.855	0.87	0.875	0.885	0.89	0.885	0.87	0.8	0.1	0.8	0.87	0.885	0.89	0.885	0.875	0.87	0.855	0.84
0.86	0.87	0.88	0.89	0.9	0.905	0.9	0.88	0.85	0.1	0.85	0.88	0.9	0.905	0.9	0.89	0.88	0.87	0.86
0.88	0.885	0.9	0.91	0.915	0.92	0.91	0.885	0.85	0.1	0.85	0.885	0.91	0.92	0.915	0.91	0.9	0.885	0.88
0.895	0.9	0.91	0.92	0.925	0.925	0.92	0.9	0.86	0.1	0.86	0.9	0.92	0.925	0.925	0.92	0.91	0.9	0.895
0.91	0.92	0.925	0.93	0.93	0.93	0.92	0.905	0.875	0.1	0.875	0.905	0.92	0.93	0.93	0.93	0.925	0.92	0.91];
								

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LIMITS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nm, maximum continuous torque corresponding to speeds in mc_map_spd
mc_max_trq=[18.36 18.36 18.15 18.04 17.94 17.83 17.73 17.62 17.57 17.17 0];
mc_max_gen_trq=-1*[18.36 18.36 18.15 18.04 17.94 17.83 17.73 17.62 17.57 17.17 0]; % estimate


mc_max_pwr_kW = max(mc_max_trq.*mc_map_spd)/1000;

mc_max_pwr_initial_kW = mc_max_pwr_kW;

% % max. overtorque (beyond continuous, intermittent operation only)
% mc_overtrq_factor=225/163;
% 
% mc_max_crrnt=300; % (A), maximum current allowed by the controller and motor
% mc_min_volts=30; % (V), minimum voltage allowed by the controller and motor
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % DEFAULT SCALING
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % (--), used to scale mc_map_spd to simulate a faster or slower running motor 
% mc_spd_scale=1.0;
% 
% % (--), used to scale mc_map_trq to simulate a higher or lower torque motor
% mc_trq_scale=1.0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OTHER DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mc_inertia=0.014;  % (kg*m^2), rotor inertia
mc_mass=14.1;  % (kg), mass of motor and controller																			

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