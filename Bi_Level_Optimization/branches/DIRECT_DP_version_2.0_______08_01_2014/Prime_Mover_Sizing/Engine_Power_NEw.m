clear all
clc
close all

% The power and energy requirement from the powertrain is determined from a given set of vehicle cruising and acceleration specifications

lb2kg = 0.453592; 
mph_mps = 0.44704;

%Vehicle Parameters
Sw = 300*lb2kg;
luggage = 30*4; % 30 kg per passenger
people = 70*4; % 70 kg per passenger
     m = 1446.96 + Sw + luggage + people ;             % kg - Total Vehicle Mass - Bass mass of a camery
    
     m =     1.9190e+03;
%      Iw = 2*4;% kg*m^2
     Jw = 4*0.19; % Inertia of the Wheel
     g = 9.81;            % m/s^2  - Gravity
     Frr = 0.01;          %Rolling Resistance
     nt = 0.9;           %Tranmission Efficiency
     rho = 1.205;            % kg/m^3  - Air Density
     Cd = 0.28;           %Aerodynamic Drag Coefficient
     Af = 2.52;            % m^2 - Frontal Area
     rwh = 0.287;           % m - Wheel Radius
     nem = 0.95;           % Motor Average Efficiency
     delta = 1+ Jw/(m*rwh^2);       % Mass Factor

     % Paux???
 % To determine the base power requirement for steady state driving

    %Base Power = For cruising of 60 mph on a flat road
        grade = 0; 
        V = 60*mph_mps;              
        Pb1 = (m*g*Frr + 0.5*rho*Cd*Af*V^2 + m*g*grade)*V/(1000*nt)

    %Base Power = For cruising of 80 mph on a flat road
        grade = 0; 
        V = 100*mph_mps;              
        Pb2 = (m*g*Frr + 0.5*rho*Cd*Af*V^2 + m*g*grade)*V/(1000*nt)

    %Base Power = For cruising of 55 mph on a road with a %5 grade
        grade = 0.05; 
        V = 40*mph_mps;              
        Pb3 = (m*g*Frr + 0.5*rho*Cd*Af*V^2 + m*g*grade)*V/(1000*nt)
 
% To determine the maximum power requirement for acceleration from 0 - 60 mph in 12 seconds
    
     %Determine Acceleration Requirement
         Vf = 60*mph_mps;               % m/s - Final Velocity
         Vi = 0*1.60934/3.6;                % m/s - Initial Velocity
         dt = 13;                           % s - time to accelearate
         dv_dt = (Vf - Vi)/dt;              % m/s^2 - Required Acceleration

     %The maximum total electric power requirement 
         Ptot = (m*g*Frr + 0.5*rho*Cd*Af*Vf^2 + m*delta*dv_dt)*Vf/(1000*nt*nem)
 
         
 % To Plot the Tractive Effort Verse the Vehicle Speed
 cd('Components');
 TEST_Run;
 Engine_41_kW;
 Motor_75_kW;
%  Vehicle_Parameters;
cd ..

V = 0:1:120*mph_mps;  
Wheel_Spd = V/rwh;   % rad/sec

%Base Power 1
grade = 0;
Peng1 = (m*g*Frr*ones(size(V)) + 0.5*rho*Cd*Af*V.^2 + m*g*grade*ones(size(V))).*V/(1000*nt);

%Base Power 2
grade = 0.05;
Peng2 = (m*g*Frr*ones(size(V)) + 0.5*rho*Cd*Af*V.^2 + m*g*grade*ones(size(V))).*V/(1000*nt);

%Base Power Z
i = 1;
for alpha = 0:5:20
F_RL(i,:) = m*g*Frr*cos(alpha*pi/180)*ones(size(V)) + 0.5*rho*Cd*Af*V.^2 + m*g*sin(alpha*pi/180)*ones(size(V));
i = i + 1;
end

x3_grid = [4.484 2.872 1.842 1.414 1.000 0.742];  % [1st 2nd...]             % Gear Level
x3_length = length(x3_grid);

C = [[1 1 0];[1 0 1];[0 1 1];[1 0 0];[0 1 0];[0 0 1]]; % Cell array of colros.
figure(1);clf
figure(2);clf
for x3 = 1:x3_length
    Gear = x3_grid(x3);
    Eng_Spd = Wheel_Spd*FD*Gear;
    Eng_Torque =  interp1(eng_map_spd,eng_max_trq,Eng_Spd);
    Eng_Power(x3,:) = Eng_Torque.*Eng_Spd;
        
    Eng_Wheel_Torque = Eng_Torque*FD*Gear*nt;  %N*m
    Eng_Tractive_Effort(x3,:) = Eng_Wheel_Torque/rwh;  % N

    figure(1);
    plot(V/mph_mps,Eng_Tractive_Effort(x3,:)/1000,'color',C(x3,:),'linewidth',5)
    hold on
    
    figure(2);
    plot(V/mph_mps,Eng_Power(x3,:)/1000,'color',C(x3,:),'linewidth',3)
    hold on

end

%Engine
Max_Eng_Tractive_Effort = max(Eng_Tractive_Effort,[],1);
figure(1);
plot(V/mph_mps,Max_Eng_Tractive_Effort/1000,'k*','markersize',8)
hold on

% Motor
Motor_Spd = Wheel_Spd*FD*G;
Motor_Torque = interp1(m_map_spd,m_max_trq,Motor_Spd);

Motor_Wheel_Torque = Motor_Torque*FD*G;
Motor_Tractive_Effort = Motor_Wheel_Torque/rwh;
plot(V/mph_mps,Motor_Tractive_Effort/1000,'k--','linewidth',3)
hold on

% Max
Max_Tractive_Effort = Max_Eng_Tractive_Effort + Motor_Tractive_Effort;
plot(V/mph_mps,Max_Tractive_Effort/1000,'k-*','linewidth',3)
hold on

for i = 1:5
figure(1);
plot(V/mph_mps,F_RL(i,:)/1000,'-ko',...
                'LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor',C(i,:),...
                'MarkerSize',5)
hold on
end
set(gca,'FontSize',20,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',25,'fontWeight','bold')
legend('Gear 1', 'Gear 2', 'Gear 3','Gear 4', 'Gear 5', 'Gear 6','Engine Max','Motor Max','Hybrid');
xlabel('Vehicle Speed (MPH)');
ylabel('Tractive Effort (kN)');
title('Road Load for Grades of: 0:5:20 degrees')
grid; hold off

figure(2);
plot(V/mph_mps,Peng1,'k','linewidth',3)
hold on
plot(V/mph_mps,Peng2,'k--','linewidth',3)
set(gca,'FontSize',20,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',25,'fontWeight','bold')
legend('Gear 1', 'Gear 2', 'Gear 3','Gear 4', 'Gear 5', 'Gear 6','Flat Road','5% Road Grade') 
xlabel('Vehicle Speed (MPH)');
ylabel('Engine Power (kW)');
grid; hold off

% Determine the Average Power Needed for the Drive Cycle
 cd('Drive_Cycle');
load CYC_HWFET; cyc_name = 'HWFET';
% load CYC_UDDS; cyc_name = 'UDDS';
% load CYC_US06; cyc_name = 'US06';
%load CYC_LA92; cyc_name = 'LA92';
%load CONST_65; cyc_name = 'CONST_65';
%load CONST_45; cyc_name = 'CONST_45';
%load CYC_COMMUTER; cyc_name = 'COMMUTER';
Manipulate_Drive_Cycle;
cd ..
%%

Pd_no_regen(Pd > 0) = Pd(Pd > 0);
Pd_ave = sum(Pd_no_regen)/length(Pd_no_regen)*ones(size(Pd));

figure(3);clf
plot(1:time_cyc,Pd/1000,'k','LineWidth',3)
hold on
plot(1:time_cyc,Pd_ave/1000,'b','LineWidth',3);
legend('Instantaneous Power','Average Power')
title(['Average Power with Zero',10,'Regenerative Braking = ',num2str(round(Pd_ave(1)/1000)),' (kW)',10,cyc_name],'fontWeight','bold','fontSize',15);
ylabel('Power (kW)')
xlabel('time (sec)','fontWeight','bold','fontSize',12);
set(gca,'fontSize',12,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',25,'fontWeight','bold'),grid

%% Vehicle speed, engine power, and resistance power vs. acceleration time
%Find the maximum engine power for every vehicle speed (From fig 2)
Max_Eng_Power = max(Eng_Power,[],1);
dm = 1.04;  % Check this!
Vf = 60*mph_mps;
Vb = 20*mph_mps;
ta = 12;
Pm = dm*m/(2*nem*ta)*(Vf^2 + Vb^2)
figure(4);
plot(V/mph_mps,Max_Eng_Power/1000,'k')
% dv_dt = Motor_Tractive_Effort*FD*G*nem/(rwh*dm*m);
% t_ = V./dv_dt
% V_ = dv_dt.*t_
% figure(5);clf
% plot(t_,V_/mph_mps)