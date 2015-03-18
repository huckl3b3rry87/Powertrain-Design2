% Define Vehicle Parameters 

% Define Load
lb2kg = 0.453592;
% Sw = 0*lb2kg;
% luggage = 10*4; % 30 kg per passenger
% people = 70*4; % 70 kg per passenger
% Load = Sw + luggage + people;

%% FROM ADVISOR %%
% Note on vehicle mass:
%		The actual average vehicle mass of a 1994 Saturn SL1 5spd is 2325 pounds.
%		If you wish to accurately set your totalvehicle mass
% 	 	to this value in the A2 GUI, you should use the mass override
%		checkbox and enter in the value 1191, which is (2325+300)/2.205 = 1191 kg, which comes from
%		adding on 300 lbs of EPA test mass, and then converting pounds to kilograms.
%		The glider mass below is just an estimate that gives 2325 pounds for a 95 kW
%		scaled SI95 engine in a conventional vehicle with 5-speed transmission.
glider_mass=(2325/2.205)-462; % (kg), vehicle mass w/o propulsion system (fuel converter,
                     % exhaust aftertreatment, drivetrain, motor, ESS, generator)
Load=136; %kg  cargo mass
gb_mass=141/2.205; % (kg), mass of the gearbox - 1990 Taurus, OTA Report
fd_mass=110/2.205; % (kg), mass of the final drive - 1990 Taurus, OTA report
tx_mass=gb_mass+fd_mass;% (kg), mass of the gearbox + final drive=(transmission)
ex_mass = 11;  % (kg)
% Base Vehicle
Base_Vehicle = glider_mass + tx_mass + ex_mass;

g = 9.81;  % (m/s^2)
rho = 1.2;     % density of air [kg/m^3]
Cd = 0.335; % Drag Coefficient

rwh = 0.33415;   % Define the radius of the tire (m)
Frr = 0.015;  % Rolling Resistance
grade = 0; % Road Grade
Af = 2;  % [m^2], vehicle frontal area (m^2) 
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