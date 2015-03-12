% Define Vehicle Parameters 

% Define Load
lb2kg = 0.453592;
Sw = 0*lb2kg;
luggage = 10*1; % 10 kg per passenger
people = 70*1; % 70 kg per passenger
Load = Sw + luggage + people;

% Base Vehicle
% Weight without Batteries	822 lb (373 kg)
% Curb Weight	1150 lb (522 kg)
% Vehicle Load Capacity	800 lb (363 kg)
Base_Vehicle = 350;

g = 9.81;  % (m/s^2)
rho = 1.2;     % density of air [kg/m^3]
Cd = 0.15; % Drag Coefficient

%http://www.etrailer.com/Tires-and-Wheels/Kenda/AM90016.html
%http://www.ezgo.com/2five/learn/specs2.php
% rwh = 5; % inches
rwh = 0.127;   % Define the radius of the tire (m)
Frr = 0.007;  % Rolling Resistance
grade = 0; % Road Grade

% Overall Width	47.3 in (120 cm)
%Overall Height	71.4 in (181 cm
Af = 1.4;  % [m^2], vehicle frontal area (m^2) 
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


Paux =  400; 

gear = [4.484 1.842 0.742];  % [1st 2nd...] 