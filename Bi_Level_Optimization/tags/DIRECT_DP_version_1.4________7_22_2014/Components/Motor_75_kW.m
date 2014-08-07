% ADVISOR Data file:  MC_AC75.m
%
% Data source:  
% Lester, L.E., et al., "An Induction Motor Power Train for EVs--The Right
% Power at the Right Price," reprinted in Proceedings:  Advanced Components for
% Electric and Hybrid Electric Vehicles, 10/27-28/93, Gaithersburg, MD.  
% Mr. Lester was employed with Westinghouse in Maryland at that time, and may
% be still.
%
% Data confidence level:  Good: data from a published paper
%
% Notes:  This is a Westinghouse, 75 kW, AC Induction motor 
%         Efficiency/loss data appropriate for a 320 V system.
%
% Created on:  6/15/98  
% By:  MRC & KW
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPEED & TORQUE RANGES over which data is defined
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (rad/s), speed range of the motor
mc_map_spd=[0 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000]*(2*pi/60);
% Note: the above conversion from RPM to rad/s was fixed 6/16/98

% (N*m), torque range of the motor
mc_map_trq=[-200 -180 -160 -140 -120 -100 -80 -60 -40 -20 ...
	0 20 40 60 80 100 120 140 160 180 200]*4.448/3.281;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EFFICIENCY AND INPUT POWER MAPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (--), efficiency map indexed vertically by mc_map_spd and 
% horizontally by mc_map_trq
mc_eff_map=[...
0.7	0.7	0.7	0.7	0.7	0.7	0.7	0.7	0.7	0.7	0.7	0.7	0.7	0.7	0.7	0.7	0.7	0.7	0.7	0.7	0.7
0.78	0.78	0.79	0.8	0.81	0.82	0.82	0.82	0.81	0.77	0.7	0.77	0.81	0.82	0.82	0.82	0.81	0.8	0.79	0.78	0.78
0.85	0.86	0.86	0.86	0.87	0.88	0.87	0.86	0.85	0.82	0.7	0.82	0.85	0.86	0.87	0.88	0.87	0.86	0.86	0.86	0.85
0.86	0.87	0.88	0.89	0.9	0.9	0.9	0.9	0.89	0.87	0.7	0.87	0.89	0.9	0.9	0.9	0.9	0.89	0.88	0.87	0.86
0.81	0.82	0.85	0.87	0.88	0.9	0.91	0.91	0.91	0.88	0.7	0.88	0.91	0.91	0.91	0.9	0.88	0.87	0.85	0.82	0.81
0.82	0.82	0.82	0.82	0.85	0.87	0.9	0.91	0.91	0.89	0.7	0.89	0.91	0.91	0.9	0.87	0.85	0.82	0.82	0.82	0.82
0.79	0.79	0.79	0.78	0.79	0.82	0.86	0.9	0.91	0.9	0.7	0.9	0.91	0.9	0.86	0.82	0.79	0.78	0.79	0.79	0.79
0.78	0.78	0.78	0.78	0.78	0.78	0.8	0.88	0.91	0.91	0.7	0.91	0.91	0.88	0.8	0.78	0.78	0.78	0.78	0.78	0.78
0.78	0.78	0.78	0.78	0.78	0.78	0.78	0.8	0.9	0.92	0.7	0.92	0.9	0.8	0.78	0.78	0.78	0.78	0.78	0.78	0.78
0.78	0.78	0.78	0.78	0.78	0.78	0.78	0.78	0.88	0.92	0.7	0.92	0.88	0.78	0.78	0.78	0.78	0.78	0.78	0.78	0.78
0.78	0.78	0.78	0.78	0.78	0.78	0.78	0.78	0.8	0.92	0.7	0.92	0.8	0.78	0.78	0.78	0.78	0.78	0.78	0.78	0.78];

%if ~exist('mc_inpwr_map')
 %  disp('Converting: MC_AC75 motor map efficiency data --> power loss data')
	%% find indices of well-defined efficiencies (where speed and torque > 0)
	pos_trqs=find(mc_map_trq>0);
	pos_spds=find(mc_map_spd>0);

	%% compute losses in well-defined efficiency area
	[T1,w1]=meshgrid(mc_map_trq(pos_trqs),mc_map_spd(pos_spds));
	mc_outpwr1_map=T1.*w1;
	mc_losspwr_map=(1./mc_eff_map(pos_spds,pos_trqs)-1).*mc_outpwr1_map;

	%% to compute losses in entire operating range
	%% ASSUME that losses are symmetric about zero-torque axis, and
     %% ASSUME that losses at zero torque are the same as those at the lowest
     %% positive torque, and
     %% ASSUME that losses at zero speed are the same as those at the lowest
     %% positive speed
	mc_losspwr_map=[fliplr(mc_losspwr_map) mc_losspwr_map(:,1) mc_losspwr_map];
	mc_losspwr_map=[mc_losspwr_map(1,:);mc_losspwr_map];

	%% compute input power (power req'd at electrical side of motor/inverter set)
	[T,w]=meshgrid(mc_map_trq,mc_map_spd);
	mc_outpwr_map=T.*w;
	mc_inpwr_map=mc_outpwr_map+mc_losspwr_map;
%end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LIMITS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% max torque curve of the motor indexed by mc_map_spd
mc_max_trq=[200 200 200 175.2 131.4 105.1 87.6 75.1 65.7 58.4 52.5]*...
	4.448/3.281; % (N*m)

mc_max_gen_trq=-1*[200 200 200 175.2 131.4 105.1 87.6 75.1 65.7 58.4 52.5]*...
	4.448/3.281; % (N*m), estimate

% maximum overtorque capability (not continuous, because the motor would overheat) 
mc_overtrq_factor=1.8; % (--), estimated

mc_max_crrnt=480; % (A), maximum current allowed by the controller and motor
mc_min_volts=120; % (V), minimum voltage allowed by the controller and motor


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
% mc_inertia=0;  % (kg*m^2), rotor inertia; unknown
% mc_mass=91;  % (kg), mass of motor and controller																			
% 
% % user definable mass scaling function
% mc_mass_scale_fun=inline('(x(1)*mc_trq_scale+x(2))*(x(3)*mc_spd_scale+x(4))*mc_mass','x','mc_spd_scale','mc_trq_scale','mc_mass');
% mc_mass_scale_coef=[1 0 1 0]; % mass scaling function coefficients

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Define all Variables for DP
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
m_map_spd = mc_map_spd;
m_map_trq = mc_map_trq;

m_max_trq = mc_max_trq;

m_eff_map = mc_eff_map;

Wm_min = -max(mc_map_spd);
Wm_max =  max(mc_map_spd);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Define all Variables for DP
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%


