% Define Vehicle Parameters 

% Define Load
lb2kg = 0.453592;
Sw = 0*lb2kg;
luggage = 10*1; % 30 kg per passenger
people = 70*1; % 70 kg per passenger
Load = Sw + luggage + people;

% Base Vehicle
Base_Vehicle = 970;

g = 9.81;  % (m/s^2)
rho = 1.2;     % density of air [kg/m^3]
Cd = 0.28; % Drag Coefficient

rwh = 0.33415;   % Define the radius of the tire (m)
Frr = 0.012;  % Rolling Resistance
grade = 0; % Road Grade
Af = 2.52;  % [m^2], vehicle frontal area (m^2) 
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


Paux =  0; 

gear = [4.484 2.872 1.842 1.414 1.000 0.742];  % [1st 2nd...]             % Gear Level