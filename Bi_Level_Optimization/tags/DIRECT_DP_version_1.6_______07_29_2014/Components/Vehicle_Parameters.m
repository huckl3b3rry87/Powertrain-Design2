% Define Vehicle Parameters 
m = 1450;  % Weight in (kg)
g = 9.81;  % (m/s^2)
rho = 1.2;     % density of air [kg/m^3]
Cd = 0.44; % Drag Coefficient

Jwh = 0.19; % Inertia of the Wheel
rwh = 0.287;   % Define the radius of the tire (m)
Frr = 0.015;  % Rolling Resistance
alpha = 0; % Road Grade (rad)
Vo =  0;  % Wind Speed (m/s)
% FD = 3.4; % Final Drive Ratio
Af = 3.87;  % [m^2], vehicle frontal area (m^2)
% G = 3;  

% Define Conversion Factors
mph_mps = 1/2.237;
rpm2rads = pi/30;
rads2rpm = 1/rpm2rads;
gasoline_density = 0.7197; % [kg/liter]
liter2gallon = 0.264172;

% Define Some battery stuff
MIN_SOC = 0.4;
MAX_SOC = 0.8;

