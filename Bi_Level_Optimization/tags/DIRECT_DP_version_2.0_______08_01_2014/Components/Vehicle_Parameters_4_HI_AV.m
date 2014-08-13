% Define Vehicle Parameters 

% Define Load
lb2kg = 0.453592;
Sw = 300*lb2kg;
luggage = 30*4; % 30 kg per passenger
people = 70*4; % 70 kg per passenger
Load = Sw + luggage + people;

% Base Vehicle
Base_Vehicle = 1121;

% Total Mass
m = Base_Vehicle + Load + fc_mass + ess_mass + mc_mass;             % kg - Total Vehicle Mass - Bass mass of a camery

g = 9.81;  % (m/s^2)
rho = 1.2;     % density of air [kg/m^3]
Cd = 0.28; % Drag Coefficient

rwh = 0.287;   % Define the radius of the tire (m)
Frr = 0.015;  % Rolling Resistance
grade = 0; % Road Grade
% FD = 3.4; % Final Drive Ratio
Af = 2.52;  % [m^2], vehicle frontal area (m^2)
% G = 3;  
nt = 1;  % Transmission Efficency

% Define Conversion Factors
mph_mps = 1/2.237;
rpm2rads = pi/30;
rads2rpm = 1/rpm2rads;
gasoline_density = 0.7197; % [kg/liter]
liter2gallon = 0.264172;

% Define Some battery stuff
MIN_SOC = 0.4;
MAX_SOC = 0.8;


Paux = 2000; 