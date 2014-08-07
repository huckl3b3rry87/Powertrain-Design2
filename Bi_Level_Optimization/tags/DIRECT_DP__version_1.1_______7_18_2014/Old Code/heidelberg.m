clear ALL
clear
close all
clc

%%
% Define Vehicle Parameters -  Eventually Define Outside of this Program
m = 1450;  % Weight in (kg)
g = 9.81;  % (m/s^2)
rho = 1.2;     % density of air [kg/m^3]
Cd = 0.28; % Drag Coefiecient
Jwh = 0.19; % Inertia of the Wheel
rwh = 0.287;   % Define the radius of the tire (m)
Frr = 0.015;  % Rolling Resistance
alpha = 0; % Road Grade (rad)
Vo =  0;  % Wind Speed (m/s)
FD = 3; % Final Drive Ratio
Af = 2.52;  % [m^2], vehicle frontal area (m^2)
G = 3;

% Define Conversion Factors
mph_mps = 1/2.237;
rpm2rads = pi/30;
rads2rpm = 1/rpm2rads;
gasoline_density = 0.7197; % [kg/liter]
liter2gallon = 0.264172;

% Define Some battery stuff
MIN_SOC = 0.4;
MAX_SOC = 0.7;

% Need to load the engine and motor
cd('Components');
Engine_2rz_0410; % with extra row for 0 torque
Motor_int;
Battery_int;

cd ..  % Come up one folder

% Define some Motor Stuff
Wm_min = -12000*rpm2rads;
Wm_max =  12000*rpm2rads;

% Define some Engine Stuff
W_eng_min = 0;
W_eng_max = max(eng_consum_spd); 

dt = 1;
load CYC_HWFET; cyc_name = 'HWFET';
% load CYC_UDDS; cyc_name = 'UDDS';
%load CYC_US06; cyc_name = 'US06';
v = cyc_mph(:,2)*mph_mps;  % Cycle Speed in (m/s)
time_cyc = length(v);
cyc_time = 1:time_cyc;
% Need to define the a(t) of the profile
for i = 2:size(v)   % Backwards difference omits first point
    a(i) = (v(i) - v(i-1))/1;
end

for t = 1:1:time_cyc  

    %% Calculate Tw, Ww and Pd from the Drive Cycle
    Tw(t) = rwh*((m*g*sin(alpha) + Frr*m*g*cos(alpha) + 0.5*rho*Cd*(v(t)+Vo)^2)*Af + (m + Jwh/rwh^2)*a(t));  % Torque needed at the wheels
    Ww(t) = (v(t)/rwh);       % Wheel Speed Needed - [rad/sec]
    Pd(t) = Ww(t)*Tw(t);      % Power Demand       - [W]
end
%%
clf
figure(1);
[AX,H1,H2] = plotyy(1:time_cyc,Tw,1:time_cyc,Ww);
set(get(AX(1),'Ylabel'),'String','Wheel Torque (Nm)','fontWeight','bold','fontSize',12)
set(get(AX(2),'Ylabel'),'String','Wheel Speed (rad/sec)','fontWeight','bold','fontSize',12)
set(H1,'LineWidth',2);
set(H2,'LineWidth',2);
set(AX(2),'fontWeight','bold','fontSize',12)
xlabel('time (sec)','fontWeight','bold','fontSize',12);
set(gca,'fontSize',12,'fontWeight','bold'),grid
title('HWFET cycle','fontWeight','bold','fontSize',16)

figure(2);
plot(1:time_cyc,cyc_mph(:,2),'LineWidth',2)
ylabel('Speed (mph)','fontWeight','bold','fontSize',12)
xlabel('time (sec)','fontWeight','bold','fontSize',12);
set(gca,'fontSize',12,'fontWeight','bold'),grid
title('HWFET cycle','fontWeight','bold','fontSize',16)

%%
%% Dynamic Programming

%% Define State Grid
x1_grid = 0.4:0.005:0.8;       % SOC 
x1_length = length(x1_grid);

%% Define Control Grids
u1_grid = [0:5:140]';          % Engine Torque [N-m]
u1_length = length(u1_grid);  

u2_grid = [4.484 2.872 1.842 1.414 1.000 0.742];  % [1st 2nd...]             % Gear Level
u2_length = length(u2_grid);

% Define Some soft cnstraints stuff
SOC_penalty = linspace(0.01,10,20);
NEAR_SOC_min = MIN_SOC + fliplr(linspace(0.001,0.02,20));

% Make some folders
tables = [cyc_name, ' TABLES'];
mkdir(tables)
cd(tables); % path into the folder

disp('-------------------------------------------------');
disp('Simulating Dynamics');
disp('-------------------------------------------------');

tic;
for t = 1:time_cyc
    table_x1 = single(zeros(x1_length,u1_length,u2_length)); % [x1]x[u1]x[u2]
    table_L = single(zeros(x1_length,u1_length,u2_length));  % [x1]x[u1]x[u2]
    
    for x1 = 1:x1_length
        SOC_c = x1_grid(x1);                  % Current SOC
        Te_c = u1_grid;                       % Engine Control, [u1]x[1] 

        for u2 = 1:u2_length
            u2_c = u2_grid(u2);               % Current Gear Level

            We_c = Ww(t)*FD*u2_c;             % [rad/sec]
            Wm_c = Ww(t)*FD*u2_c*G;             % [rad/sec]

            if Pd(t) > 0 % Driving
                Tm_c = Tw(t)/(u2_c*FD*G)*ones(size(Te_c)) - Te_c;  % [u1]x[1]

            else         % Braking
                Te_c = zeros(u1_length,1);   % Turn The Engine Off
                We_c = 0;
                Tm_c = Tw(t)/(u2_c*FD*G)*ones(size(Te_c)) - Te_c;  % [u1]x[1]
            end
             
            % Check Motor
            Tm_max(:,u2) = interp1(m_map_spd,m_max_trq,abs(Wm_c))*ones(size(Tm_c));       
            Tm_save(:,u2) = Tm_c;                                                   
            Wm_save(:,u2) = Wm_c*ones(size(u1_grid));                    %  Check For Each Gear
            
            % Check Engine
            We_save(:,u2) = We_c*ones(size(u1_grid));                    % Check For Each Gear
            
            % Update x1 
            eff_m = interp2(m_map_trq, m_map_spd, m_eff_map, Tm_c, abs(Wm_c))';   % [u1]x[1]  - should be this size
            eff_m(isnan(eff_m))=0.2; % to avoid nan
            Pm_c = (Wm_c*Tm_c).*(eff_m.^-sign(Wm_c*Tm_c));
            Pbat = Pm_c;     % [u1]x[1]  - Need to eventually check Pbat
            
            if Pbat>=0
                Pbatt_max = interp1(ess_soc, ess_max_pwr_dis, SOC_c);
                Pbat(Pbat>Pbatt_max) = Pbatt_max;
                rint_c = interp1(ess_soc,ess_r_dis,SOC_c);
            else
                Pbatt_min = interp1(ess_soc, ess_max_pwr_chg, SOC_c);
                Pbat(Pbat<-Pbatt_min) = -Pbatt_min;
                rint_c = interp1(ess_soc,ess_r_chg,SOC_c);
            end
            voc_c = interp1(ess_soc,ess_voc,SOC_c);
            SOC_n = real((-(voc_c-(voc_c.^2-4*Pbat.*rint_c).^(1/2))./(2*rint_c)/(ess_cap_ah*3600))*dt + SOC_c); % [u1]x[1]
            fuel = interp2(eng_consum_trq,eng_consum_spd,eng_consum_fuel,Te_c,We_c,'linear')*dt;   %[u1]x[1]
     
            SOC_save(:,u2) = SOC_n;                       %[u1]x[u2] - Nine Different SOC's depending on the control
            inst_fuel(:,u2) = fuel;                       %[u1]x[u2] - Nine Different fuel's depending on the control
            
        end  % End of u2 ( Gear Level ) Loops
        
        % Check Motor     
        infeasible_Tm = (Tm_save > Tm_max) | (Tm_save < -Tm_max);        %[u1]x[u2]      
        infeasible_Wm = (Wm_save > Wm_max) | (Wm_save < Wm_min);
        
        % Check Engine - Torque is defined, so it will be ok. If we had engine speed as a state, we could check the acceleration
        infeasible_We = (We_save < 0) | (We_save > W_eng_max);                 

        % Check SOC
        infeasible_SOC = (SOC_save < MIN_SOC) | (SOC_save > MAX_SOC);  % [u1]x[u2]
        SOC_save(SOC_save > MAX_SOC) = MAX_SOC;  % Why saturate, isn't that cheating ?
        SOC_save(SOC_save < MIN_SOC) = MIN_SOC;

        table_x1(x1,:,:) = SOC_save;   % Next state (SOC), for this particualar SOC grid point and time step,  [x1]x[u1]x[u2]
        table_L(x1,:,:) = inst_fuel; 
        
        % Constraints
        table_L(x1,infeasible_We) = 10;
        table_L(x1,infeasible_Tm) = 10;
        table_L(x1,infeasible_Wm) = 10; 
        table_L(x1,infeasible_SOC) = 10;

        SOC_soft_constraint = zeros(size(SOC_save));
        SOC_soft_constraint(NEAR_SOC_min(1)>SOC_n) =  SOC_penalty(1);
        SOC_soft_constraint(NEAR_SOC_min(2)>SOC_n) =  SOC_penalty(2);
        SOC_soft_constraint(NEAR_SOC_min(3)>SOC_n) =  SOC_penalty(3);
        SOC_soft_constraint(NEAR_SOC_min(4)>SOC_n) =  SOC_penalty(4);
        SOC_soft_constraint(NEAR_SOC_min(5)>SOC_n) =  SOC_penalty(5);
        SOC_soft_constraint(NEAR_SOC_min(6)>SOC_n) =  SOC_penalty(6);
        SOC_soft_constraint(NEAR_SOC_min(7)>SOC_n) =  SOC_penalty(7);
        SOC_soft_constraint(NEAR_SOC_min(8)>SOC_n) =  SOC_penalty(8);
        SOC_soft_constraint(NEAR_SOC_min(9)>SOC_n) =  SOC_penalty(9);
        SOC_soft_constraint(NEAR_SOC_min(10)>SOC_n) =  SOC_penalty(10);
        SOC_soft_constraint(NEAR_SOC_min(11)>SOC_n) =  SOC_penalty(11);
        SOC_soft_constraint(NEAR_SOC_min(12)>SOC_n) =  SOC_penalty(12);
        SOC_soft_constraint(NEAR_SOC_min(13)>SOC_n) =  SOC_penalty(13);
        SOC_soft_constraint(NEAR_SOC_min(14)>SOC_n) =  SOC_penalty(14);
        SOC_soft_constraint(NEAR_SOC_min(15)>SOC_n) =  SOC_penalty(15);
        SOC_soft_constraint(NEAR_SOC_min(16)>SOC_n) =  SOC_penalty(16);
        SOC_soft_constraint(NEAR_SOC_min(17)>SOC_n) =  SOC_penalty(17);
        SOC_soft_constraint(NEAR_SOC_min(18)>SOC_n) =  SOC_penalty(18);
        SOC_soft_constraint(NEAR_SOC_min(19)>SOC_n) =  SOC_penalty(19);
        SOC_soft_constraint(NEAR_SOC_min(20)>SOC_n) =  SOC_penalty(20);

        temp = squeeze(table_L(x1,:,:));
        table_L(x1,:,:) = temp + SOC_soft_constraint;
        
    end  % End of x1 loop
        
    savename = ['Transitional Cost = ',num2str(t),' Table.mat'];
    save(savename,'table_x1','table_L');
    toc
end
cd ..  % Come out of the folder           
%%        
disp('-------------------------------------------------');
disp('Dynamic Programming');
disp('-------------------------------------------------');

tic;
clear table_*

% Define Parameters
BETA = 4000;
Desired_SOC = 0.55; 

folder = [cyc_name, ' TABLES'];
cd(folder);                     % Go into the tables that were previously made

%
J_STAR = repmat(BETA*(x1_grid - Desired_SOC).^2, 1, 1); % terminal state penalty, [x1]x[1]
for t = time_cyc:-1:1   
    loadfile_name = ['Transitional Cost = ',num2str(t),' TABLE.mat'];
    load(loadfile_name);
    
    J_temp = table_L + interp1(x1_grid,J_STAR,table_x1,'linear'); % [x1]x[u1]x[u2] - u2 is tileing - interp is the transitional cost at the new interpolated state grid
    [opt_value_u1, opt_id_u1] = min(J_temp,[],2);  % You are trying to find the optimal control for all states!!! This is why it is along this dimension

    [opt_value,opt_id_u2] = min(opt_value_u1,[],3);  % Check Dimensions
    opt_id_u1 = squeeze(opt_id_u1);
    opt_id_u1 =  opt_id_u1(:,1);  % Check this 
    opt_id_u2 = squeeze(opt_id_u2);
    
    J_STAR = opt_value;   % Next terminal cost!
    if all(isnan(J_STAR(:)))
        disp('NO Solution!!')
    end
    
    savename=['Cost & Control = ',num2str(t),' TABLE.mat'];
    save(savename,'J_STAR','opt_id_u1','opt_id_u2');
end
toc;
cd ..

disp('-------------------------------------------------');
disp('Simulating Final Run');
disp('-------------------------------------------------');

folder = [cyc_name, ' TABLES'];
cd(folder);

% Pre allocate for speed
SOC_final = zeros(1, time_cyc);
inst_fuel = zeros(1,time_cyc);
W_mot = zeros(1, time_cyc);
T_mot = zeros(1, time_cyc);
W_eng = zeros(1, time_cyc);
T_eng = zeros(1,time_cyc);

SOC_c = 0.55;  % Initialize SOC

tic;
for t = 1:1:time_cyc
    
    load(['Cost & Control = ',num2str(t),' TABLE.mat']);
    id_lookup_u2 = interp1(x1_grid, opt_id_u2, SOC_c,'nearest');
    
    if isnan(id_lookup_u2)
        if SOC_c < MIN_SOC
            id_lookup_u2 = opt_id_u2(1); % Check later to make sure this corresponds with min SOC
        else
            id_lookup_u2 = opt_id_u2(u2_length);
        end
    end
            
    u2_c = u2_grid(id_lookup_u2);               % Current Gear Level - check the min thing - should all be the same
     
    We_c = Ww(t)*FD*u2_c;             % [rad/sec]
    Wm_c = Ww(t)*FD*u2_c*G;           % [rad/sec]

    if Pd(t) <= 0 % Stopped or Braking
        Te_c = 0;   % Turn The Engine Off
        We_c = 0;
        Tm_c = Tw(t)/(u2_c*FD*G) - Te_c;  % [1]x[1]
        
    else         % Driving
        id_lookup_u1 = interp1(x1_grid, opt_id_u1, SOC_c,'nearest');
        
        if isnan(id_lookup_u1)
            if SOC_c < MIN_SOC
                id_lookup_u1 = opt_id_u1(1); % Check later to make sure this corresponds with min SOC
            else
                id_lookup_u1 = opt_id_u1(u1_length);
            end
        end
        Te_c = u1_grid(id_lookup_u1);  % Just grab a row of u1 eventulally
        Tm_c = Tw(t)/(u2_c*FD*G) - Te_c;  % [1]x[1]
    end  
    Tm_max = interp1(m_map_spd,m_max_trq,abs(Wm_c)); % [1]x[1]
    Tm_c(Tm_c>Tm_max) = Tm_max;
    Tm_c(Tm_c<-Tm_max) = -Tm_max;
    
    % Update x1 
    eff_m = interp2(m_map_trq, m_map_spd, m_eff_map, Tm_c, abs(Wm_c));
    eff_m(isnan(eff_m))=0.2; % to avoid nan
    Pm_c = (Wm_c.*Tm_c).*(eff_m.^-sign(Wm_c.*Tm_c));
    
    Pbat = Tm_c*Wm_c;     % [u1]x[1]
    
    if Pbat>=0
        Pbatt_max = interp1(ess_soc, ess_max_pwr_dis, SOC_c);
        Pbat(Pbat>Pbatt_max) = Pbatt_max;
        rint_c = interp1(ess_soc,ess_r_dis,SOC_c);
    else
        Pbatt_min = interp1(ess_soc, ess_max_pwr_chg, SOC_c);
        Pbat(Pbat<-Pbatt_min) = -Pbatt_min;
        rint_c = interp1(ess_soc,ess_r_chg,SOC_c);
    end
    voc_c = interp1(ess_soc,ess_voc,SOC_c);
    
    SOC_n = real((-(voc_c-(voc_c.^2-4*Pbat.*rint_c).^(1/2))./(2*rint_c)/(ess_cap_ah*3600))*dt + SOC_c);
    fuel = interp2(eng_consum_trq,eng_consum_spd,eng_consum_fuel,Te_c,We_c,'linear')*dt;   %[1]x[1]
    
    W_eng(t) = We_c;
    T_eng(t) = Te_c;
    W_mot(t) = Wm_c;
    T_mot(t) = Tm_c;
    inst_fuel(t) = fuel;
    SOC_final(t) = SOC_c;
    GEAR(t) = id_lookup_u2;
    SOC_c = SOC_n;
end
toc;

T_eng(T_eng<0) = 0;
dSOC = SOC_final(end) - SOC_final(1)

total_fuel_gram = sum(inst_fuel);
total_distance_mile = sum(cyc_mph(:,2))/3600;
mpg = total_distance_mile/(total_fuel_gram/1000/gasoline_density*liter2gallon)
%
cd ..    
%%    
% ================================================================== %%
figure(1); clf;
set(gcf, 'units', 'inch', 'pos', [7.7604    1.1146    5.8333    8.2813]);
subplot(5,1,1);
[haxes,hline1,hline2] = plotyy(cyc_time, v, cyc_time,SOC_final);
set(hline1,'marker','x', 'markersize', 3, 'markerf', 'b');
set(hline2,'marker','x', 'markersize', 3, 'markerf', [0 0.5 0]);
set(haxes, 'fontsize', 8);
set(haxes(1), 'ylim', [0 80], 'ytick', 0:20:80, 'box', 'off', 'yminortick', 'on');
set(haxes(2), 'ylim', [0.4 0.7], 'ytick', 0.4:0.1:0.7, 'box', 'off', 'yminortick', 'on');
ylabel('Speed (mph)', 'parent', haxes(1));
ylabel('SOC (-)', 'parent', haxes(2));
set(haxes(2), 'XAxisLocation','top', 'xtick', [], 'xcolor', 'w');
set(gca, 'xgrid', 'on');
title(['\alpha = ', num2str(BETA), '; SOC_{ini} = ', num2str(SOC_final(1))]);

% ==============================
subplot(5,1,2);
plot(cyc_time, inst_fuel, 'kx-', 'markersize', 3);
ylim([0 5]);
set(gca, 'pos', [0.1282    0.5935    0.7750    0.1495]);
set(gca, 'fontsize', 8);
set(gca, 'xgrid', 'on');
ylabel('Fuel Rate(g/s)');

% ==============================
subplot(5,1,3);
plot(cyc_time, W_eng*rads2rpm, 'kx-', ...
     cyc_time, W_mot*rads2rpm, 'rx-', 'markersize', 3);
set(gca, 'pos', [0.1282    0.4171    0.7750    0.1495]);
ylim([-5000 15000]);
set(gca, 'fontsize', 8);
set(gca, 'xgrid', 'on');
ylabel('Speed (rpm)');
legend('Engine', 'MG');

% ==============================
subplot(5,1,4);
plot(cyc_time, T_eng, 'kx-', ...
     cyc_time, T_mot, 'rx-', 'markersize', 3);
set(gca, 'pos', [0.1282    0.2408    0.7750    0.1495]);
ylim([-100 150]);
set(gca, 'fontsize', 8);
set(gca, 'xgrid', 'on');
ylabel('Torque (Nm)');

% ==============================
subplot(5,1,5);
plot(cyc_time, GEAR, 'kx-','markersize', 3);
set(gca, 'fontsize', 8);
set(gca, 'xgrid', 'on');
ylabel('Gear');


    
%%
figure(2);clf;
[AX,H1,H2] = plotyy(cyc_time,v,cyc_time,SOC_final);
set(get(AX(1),'Ylabel'),'String','Vehicle Speed (mph)','fontWeight','bold','fontSize',20)
set(get(AX(2),'Ylabel'),'String','SOC','fontWeight','bold','fontSize',20)
set(H1,'LineWidth',2);
set(H2,'LineWidth',2);
set(AX(2),'fontWeight','bold','fontSize',20)
xlabel('time (sec)','fontWeight','bold','fontSize',20);
set(gca,'fontSize',20,'fontWeight','bold'),grid
% title('HWFET cycle','fontWeight','bold','fontSize',20)
       
%%
figure(3);clf;
[AX,H1,H2] = plotyy(cyc_time,v,cyc_time,GEAR);
set(get(AX(1),'Ylabel'),'String','Vehicle Speed (mph)','fontWeight','bold','fontSize',20)
set(get(AX(2),'Ylabel'),'String','Optimal Gear Level','fontWeight','bold','fontSize',20)
set(H1,'LineWidth',2);
set(H2,'LineWidth',2);
set(AX(2),'fontWeight','bold','fontSize',20)
xlabel('time (sec)','fontWeight','bold','fontSize',20);
set(gca,'fontSize',20,'fontWeight','bold'),grid
% title('HWFET cycle','fontWeight','bold','fontSize',20)

%%
figure(4);clf;
[AX,H1,H2] = plotyy(cyc_time,v,cyc_time,inst_fuel);
set(get(AX(1),'Ylabel'),'String','Vehicle Speed (mph)','fontWeight','bold','fontSize',20)
set(get(AX(2),'Ylabel'),'String','Instantaneous Fuel (g/s)','fontWeight','bold','fontSize',20)
set(H1,'LineWidth',2);
set(H2,'LineWidth',2);
set(AX(2),'fontWeight','bold','fontSize',20)
xlabel('time (sec)','fontWeight','bold','fontSize',20);
set(gca,'fontSize',20,'fontWeight','bold'),grid
% title('HWFET cycle','fontWeight','bold','fontSize',20)

%%
figure(5);clf;
plot(cyc_time,W_eng*rads2rpm,'color','g','LineWidth',2);
hold on 
plot(cyc_time,W_mot*rads2rpm,'LineWidth',2);grid on;
legend('Engine Speed','Motor Speed')
ylabel('Speed (rpm)','fontWeight','bold','fontSize',20)
xlabel('time (sec)','fontWeight','bold','fontSize',20);
set(gca,'fontSize',20,'fontWeight','bold')
% title('HWFET cycle','fontWeight','bold','fontSize',20)

%%
figure(6);clf;
plot(cyc_time,T_eng,'color','g','LineWidth',2);
hold on 
plot(cyc_time,T_mot,'LineWidth',2);grid on;
legend('Engine Torque (Nm)','Motor Torque (Nm)')
ylabel('Torque (Nm)','fontWeight','bold','fontSize',20)
xlabel('time (sec)','fontWeight','bold','fontSize',20);
set(gca,'fontSize',20,'fontWeight','bold')
% title('HWFET cycle','fontWeight','bold','fontSize',20)

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
