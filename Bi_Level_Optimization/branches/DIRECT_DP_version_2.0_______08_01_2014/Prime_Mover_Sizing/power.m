clear all
clc
close all

% The power and energy requirement from the powertrain is determined from a given set of vehicle cruising and acceleration specifications

lb2kg = 0.453592; 
mph2mps = 0.44704;

%Vehicle Parameters
Sw = 300*lb2kg;
luggage = 30*4; % 30 kg per passenger
people = 70*4; % 70 kg per passenger
     m = 1446.96 + Sw + luggage + people ;             % kg - Total Vehicle Mass - Bass mass of a camery
    
     m =    1.9013e+03;
     
     g = 9.81;            % m/s^2  - Gravity
     fr = 0.015;          %Rolling Resistance
     nt = 0.92;           %Tranmission Efficiency
     pa = 1.2;            % kg/m^3  - Air Density
     Cd = 0.28;           %Aerodynamic Drag Coefficient
     Af1 = 78/12*77/12*0.3048*0.3048;
     Af = 0.33;            % m^2 - Frontal Area
     r = 0.287;           % m - Wheel Radius
     nem = 0.85;           % Motor Average Efficiency
     delta = 1.035;       % Mass Factor

     
 % To determine the base power requirement for steady state driving

    %Base Power = For cruising of 60 mph on a flat road
        grade = 0; 
        V = 60*mph2mps;               % km/h - vehicle Speed
        Pb1 = (m*g*fr + 0.5*pa*Cd*Af*V^2 + m*g*sin(grade))*V/(1000*nt*nem)

    %Base Power = For cruising of 80 mph on a flat road
        grade = 0; 
        V = 80*mph2mps;               % km/h - vehicle Speed
        Pb2 = (m*g*fr + 0.5*pa*Cd*Af*V^2 + m*g*sin(grade))*V/(1000*nt*nem)

    %Base Power = For cruising of 55 mph on a road with a %5 grade
        grade = 0.05; 
        V = 55*mph2mps;               % km/h - vehicle Speed
        Pb3 = (m*g*fr + 0.5*pa*Cd*Af*V^2 + m*g*sin(grade))*V/(1000*nt*nem)
 
% To determine the maximum power requirement for acceleration from 0 - 60 mph in 12 seconds
    
     %Determine Acceleration Requirement
         Vf = 60*mph2mps;               % m/s - Final Velocity
         Vi = 0*1.60934/3.6;                % m/s - Initial Velocity
         dt = 13;                           % s - time to accelearate
         dv_dt = (Vf - Vi)/dt;              % m/s^2 - Required Acceleration

     %The maximum total electric power requirement 
         Ptot = (m*g*fr + 0.5*pa*Cd*Af*Vf^2 + m*delta*dv_dt)*Vf/(1000*nt*nem)
 
 
    