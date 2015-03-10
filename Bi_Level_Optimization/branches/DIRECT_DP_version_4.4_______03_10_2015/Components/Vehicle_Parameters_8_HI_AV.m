% Define Vehicle Parameters 
% sources: http://www.edmunds.com/dodge/caravan/2007/features-specs.html?style=100809423
%          http://mcgrefer.com/sizeinfo/2156516
% Define Load
lb2kg = 0.453592;
Sw = 300*lb2kg;
luggage = 30*8; % 30 kg per passenger
people = 70*8; % 70 kg per passenger
Load = Sw + luggage + people;

% Base Vehicle
Base_Vehicle =   1.8652e+03;

g = 9.81;  % (m/s^2)
rho = 1.2;     % density of air [kg/m^3]
Cd = 0.35; % Drag Coefficient

rwh = 0.34295;   % Define the radius of the tire (m)
Frr = 0.015;  % Rolling Resistance
grade = 0; % Road Grade
% 6 ft 6.9in x 5 ft 8.9 in
% 78.9 in x 68.9 in =    3.5072 m^2
Af = 3.51;  % [m^2], vehicle frontal area (m^2) 
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


Paux = 2500 + 700; 

gear = [4.484 2.872 1.842 1.414 1.000 0.742];  % [1st 2nd...]             % Gear Level