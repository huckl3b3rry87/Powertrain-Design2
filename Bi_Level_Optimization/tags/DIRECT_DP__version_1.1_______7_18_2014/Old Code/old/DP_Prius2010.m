clear all;
close all;
clc;
format compact;

dt = 1;


%% ================================================================== %%
%% unit conversion:
rpm2rads = pi/30;
rads2rpm = 1/rpm2rads;
mph2mps = 1609/3600;
mps2mph = 1/mph2mps;
g2gallon=1/841.4/3.785;

gasoline_density = 0.7197; % [kg/liter]
liter2gallon = 0.264172;


%% ================================================================== %%
%% component initialization
cd('Components');
Vehicle_int;
Engine_2rz_0410; % with extra row for 0 torque
Motor_int;
Battery_int;
cd ..

%% engine
We_min = 0;
We_max = max(eng_consum_spd); % i.e. 5200 rpm

%% motor
mg1_max_spd = g_max_spd;
mg1_max_trq = g_max_trq;
mg1_map_spd = g_map_spd;
mg1_map_trq = g_map_trq;
mg1_eff_map = g_eff_map;
Wmg1_min = -12000*rpm2rads;
Wmg1_max =  12000*rpm2rads;

%% generator
mg2_max_spd = m_max_spd;
mg2_max_trq = m_max_trq;
mg2_map_spd = m_map_spd;
mg2_map_trq = m_map_trq;
mg2_eff_map = m_eff_map;
Wmg2_min = -12000*rpm2rads;
Wmg2_max =  12000*rpm2rads;

%% power
Pmg1_max = 40000;
Pmg2_max = 40000;
soc_min = 0.4;
soc_max = 0.7;


%% ================================================================== %%
%% powertrain parameters
%% planetary gearset parameters
R = 2.6;
S = 1;

FR = 3.3; % final drive ration of 2010 Prius (third generation Prius)
MR = 2.63; % fixed ratio on MG2 in 2010 Prius (mentioned in Namwook Kim's journal paper)

I_r=0; % ring gear inertia
I_s=0; % sun gear inertia
I_c=0; % carrier gear inertia

g_inertia=0.0000226;
m_inertia=0.015;

I_mg1 = g_inertia;% generator inertia
I_mg2 = m_inertia;% motor inertia
I_e = e_inertia; % engine inertia
I_v = M_total*R_tire^2/FR^2; % vehicle inertia


%% ================================================================== %%
%% Prius configuration
D = [R+S -R -S];

%% Mode 4: split mode
D_mode4 = [I_e+I_c, 0,         0,         D(1); % engine (carrier gear)
           0,       I_v+I_mg2, 0,         D(2); % wheel & mg2 (ring gear)
           0,       0,         I_mg1+I_s, D(3); % mg1 (sun gear)
           D(1),    D(2),      D(3),      0];
T2W = inv(D_mode4);


%% ================================================================== %%
%% Load driving cycle
cd('CYC');

load nedc.mat; cyc_name = 'necd';
cyc_time = 1:length(sch_cycle);
% cyc_time = 1:200; % only take the first 200 seconds
cyc_mps = sch_cycle(cyc_time,2); % [m/s]
cyc_mph = cyc_mps*mps2mph;

% load CYC_FUDS.mat; cyc_name = 'FUDS';
% cyc_time = 1:length(cyc_mph);
% cyc_mps = cyc_mph(:,2)*mph2mps;
% cyc_mph = cyc_mph(:,2);

% load CYC_US06.mat; cyc_name = 'US06';
% cyc_time = 1:length(cyc_mph);
% cyc_mps = cyc_mph(:,2)*mph2mps;
% cyc_mph = cyc_mph(:,2);

cyc_mps(cyc_mps<0) = 0;

cd ..

%% ==============================
mps = cyc_mps';
acc = [diff(mps)  0];% [m/s^2]

%% Drag forces, [N]
resis_acc  = M_total * acc;
resis_roll = M_total*g_gravity*f_rolling.*(mps>0);
resis_aero = 1/2*C_d*rho_air*A_frontal*mps.^2;

Wv = mps./R_tire; % [rad/s], wheel rotational speed
Tv = R_tire * (resis_aero + resis_roll + resis_acc); % +: driving, -:regen
Pv = Wv.*Tv;
PowerD = Pv; % +: driving, -:regen

%% Output node on the planetary gearset (i.e. ring gear for Prius)
gr_eff=1;
To = (R_tire * (resis_aero + resis_roll))/FR; % +: driving, -:regen (note that vehicle acceleration is not included!!!)
Wo = Wv*FR; % [rad/s], output gear speed
dWo = [diff(Wo), 0]; % [rad/s^s], output gear rotational acceleartion


%% ================================================================== %%
%% Dynamic Programming:
%%  Step1: calculate transitional cost (~415 seconds)
%%  Step2: calculate cumulative cost (~135 seconds)
%%  Step3: retrieve optimal solution (~5.6 seconds)

soc_des = 0.55; % desired SOC
alpha = 50000 % weighting coefficient for penalty on SOC deviation

%% 1st state: Engine speed
x1_We_grid = (0:100:5200)*rpm2rads; % [rad/s]
N_x1_We = length(x1_We_grid);

%% 2nd state: SOC
x2_SOC_grid = 0.4:0.005:0.7;
N_x2_SOC = length(x2_SOC_grid);

%% 1st control: Tmg1
u1_Tmg1_grid = -140:5:140; % [Nm]
N_u1_Tmg1 = length(u1_Tmg1_grid);


disp('===============================================================');
disp('Starging Step 1');
tic;
tablefolder_name = [cyc_name, '_table_folder'];
mkdir(tablefolder_name)
cd(tablefolder_name); % path into the folder
for t = 1:length(cyc_time)
    table_x1n = single(zeros(N_x1_We,N_x2_SOC,N_u1_Tmg1)); % [x1]x[x2]x[u1]
    table_x2n = single(zeros(N_x1_We,N_x2_SOC,N_u1_Tmg1)); % [x1]x[x2]x[u1]
    table_L = single(zeros(N_x1_We,N_x2_SOC,N_u1_Tmg1)); % [x1]x[x2]x[u1]
    
    Wo_c = Wo(t); % current output speed (rad/s)
    To_c = To(t);
    dWo_c = dWo(t);
    
    for a = 1:N_x1_We % loop for x1
        We_c = x1_We_grid(a);
        
        if PowerD(t) > 0
            Te_c = interp1(We_optbsfc, Te_optbsfc, We_c);
            Tmg1_c = u1_Tmg1_grid; % [1]x[u1], vector
        else % PowerD(t) <=0
            Te_c = -6; % engine drag, which slows down the engine speed
            Tmg1_c = zeros(1, N_u1_Tmg1);
        end
        
        %% generator
        Wmg1_c = (We_c*D(1) + Wo_c*D(2))/(-D(3)); % [1]x[1], scaler
        infeasible_Wmg1 = (Wmg1_c>Wmg1_max) | (Wmg1_c<Wmg1_min); % [1]x[1]
        %Wmg1_c(Wmg1_c>Wmg1_max) = Wmg1_max;
        %Wmg1_c(Wmg1_c<Wmg1_min) = Wmg1_min;
        
        %% motor
        Tmg2_c = (((dWo_c - T2W(2,1)*Te_c - T2W(2,3).*Tmg1_c)./T2W(2,2) + To_c))/MR; % [1]x[u1]
        Wmg2_c = Wo_c*MR; % [1]x[1]
        Tmg2_max = interp1(m_map_spd,m_max_trq,abs(Wmg2_c)); % [1]x[1]
        infeasible_Tmg2 = (Tmg2_c>Tmg2_max) | (Tmg2_c<-Tmg2_max); % [1]x[u1]
        %Tmg2_c(Tmg2_c>Tmg2_max) = Tmg2_max;
        %Tmg2_c(Tmg2_c<-Tmg2_max) = -Tmg2_max;
        
        %% update x1 (engine speed)
        We_n = (T2W(1,1)*Te_c + T2W(1,2).*(Tmg2_c*MR - To_c) + T2W(1,3).*Tmg1_c)*dt + We_c; % [1]x[u1]
        dWe_n = We_n - We_c;
        infeasible_We = (We_n<0) | (We_n>We_max) | (abs(dWe_n) > 400*rpm2rads); % [1]x[u1]
        We_n(We_n<0) = 0;
        %We_n(We_n>We_max) = We_max;
        
        fuel_c = dt*0.5*(interp2(eng_consum_trq,eng_consum_spd,eng_consum_fuel,Te_c,We_c,'linear') + ...
            interp2(eng_consum_trq,eng_consum_spd,eng_consum_fuel,Te_c,We_n,'linear')'); % [1]x[u1]
        
        %% update x2 (battery SOC)
        eff_mg1 = interp2(mg1_map_trq, mg1_map_spd, mg1_eff_map, Tmg1_c, abs(Wmg1_c));
        eff_mg2 = interp2(mg2_map_trq, mg2_map_spd, mg2_eff_map, Tmg2_c, abs(Wmg2_c));
        eff_mg1(isnan(eff_mg1))=0.2; % to avoid nan
        eff_mg2(isnan(eff_mg2))=0.2; % to avoid nan
        
        Pmg1_c = (Wmg1_c.*Tmg1_c).*(eff_mg1.^-sign(Wmg1_c.*Tmg1_c));
        Pmg2_c = (Wmg2_c.*Tmg2_c).*(eff_mg2.^-sign(Wmg2_c.*Tmg2_c));
        
        Pbatt = Pmg1_c + Pmg2_c; % +: discharge; -: charge; [1]x[u1]

        Pbatt = repmat(Pbatt, N_x2_SOC, 1); % [x2]x[u1]
        
        soc_c = repmat(x2_SOC_grid', 1, N_u1_Tmg1); % [x2]x[u1]
        voc_c = interp1(ess_soc,ess_voc,soc_c); % [x2]x[u1]
        Pbatt_max = interp1(ess_soc, ess_max_pwr_dis, soc_c);
        Pbatt(Pbatt>Pbatt_max) = Pbatt_max(Pbatt>Pbatt_max);
        Pbatt_min = interp1(ess_soc, ess_max_pwr_chg, soc_c);
        Pbatt(Pbatt<-Pbatt_min) = -Pbatt_min(Pbatt<-Pbatt_min);
        %infeasible_Pbatt = (Pbatt<-Pbattery_max) | (Pbatt>Pbattery_max); % [x2]x[u1]
        
        rint_c_dis = interp1(ess_soc,ess_r_dis,soc_c); % [x2]x[u1]
        rint_c_chg = interp1(ess_soc,ess_r_chg,soc_c);
        rint_c = rint_c_dis;
        rint_c(Pbatt<0) = rint_c_chg(Pbatt<0);
        
        %I = (voc_c-sqrt(voc_c.^2-4*Pbatt.*rint_c))./(2*rint_c); % battery discharge current
        soc_n = (-(voc_c-(voc_c.^2-4*Pbatt.*rint_c).^(1/2))./(2*rint_c)/(ess_cap_ah*3600))*dt + soc_c; % [x2]x[u1]
        infeasible_soc = (soc_n<soc_min) | (soc_n>soc_max); % [x2]x[u1]
        soc_n(soc_n>soc_max) = soc_max;
        soc_n(soc_n<soc_min) = soc_min;
        
        table_x1n(a,:,:) = repmat(We_n, N_x2_SOC, 1); % [x2]x[u1]
        table_x2n(a,:,:) = soc_n; % [x2]x[u1]
        table_L(a,:,:) = repmat(fuel_c, N_x2_SOC, 1); % [x2]x[u1]
        
        %% penalize infeasible operation with FINITE penalty
        table_L(a,:,infeasible_We) = 10;
        table_L(a,:,infeasible_Tmg2) = 10;
        if infeasible_Wmg1==1
            table_L(a,:,:) = 10;
        end
        table_L(a,infeasible_soc) = 10;
        
        soc_softcstr = zeros(size(soc_n));
        soc_softcstr(soc_n<(soc_min+0.02)) = 0.001;
        soc_softcstr(soc_n<(soc_min+0.019)) = 0.003;
        soc_softcstr(soc_n<(soc_min+0.018)) = 0.005;
        soc_softcstr(soc_n<(soc_min+0.017)) = 0.0068;
        soc_softcstr(soc_n<(soc_min+0.016)) = 0.0122;
        soc_softcstr(soc_n<(soc_min+0.015)) = 0.022;
        soc_softcstr(soc_n<(soc_min+0.014)) = 0.0338;
        soc_softcstr(soc_n<(soc_min+0.013)) = 0.05;
        soc_softcstr(soc_n<(soc_min+0.012)) = 0.068;
        soc_softcstr(soc_n<(soc_min+0.011)) = 0.1020;
        soc_softcstr(soc_n<(soc_min+0.01)) = 0.1650;
        soc_softcstr(soc_n<(soc_min+0.009)) = 0.2980;
        soc_softcstr(soc_n<(soc_min+0.008)) = 0.4600;
        soc_softcstr(soc_n<(soc_min+0.007)) = 0.7000;
        soc_softcstr(soc_n<(soc_min+0.006)) = 1.2;
        soc_softcstr(soc_n<(soc_min+0.005)) = 2.14;
        soc_softcstr(soc_n<(soc_min+0.004)) = 3.38;
        soc_softcstr(soc_n<(soc_min+0.003)) = 5.2;
        soc_softcstr(soc_n<(soc_min+0.002)) = 7.33;
        soc_softcstr(soc_n<(soc_min+0.001)) = 10;
        temp = squeeze(table_L(a,:,:));
        table_L(a,:,:) = temp + soc_softcstr;

    end % end of loop a
    
    savename = ['TRSN_t=',num2str(t),'_Table.mat'];
    save(savename,'table_x1n','table_x2n','table_L');
    toc
end
cd ..


disp('===============================================================');
disp('Starging Step 2');
tic;
clear table_*

target_folder = [cyc_name, '_table_folder'];
cd(target_folder);

J_star = repmat(alpha*(x2_SOC_grid - soc_des).^2, N_x1_We, 1); % terminal state penalty, [x1]x[x2]
for t = length(cyc_time):-1:1   
    loadfile_name = ['TRSN_t=',num2str(t),'_Table.mat'];
    load(loadfile_name);
    
    J_temp = table_L + interp2(x2_SOC_grid,x1_We_grid,J_star,table_x2n,table_x1n,'linear'); % [x1]x[x2]x[u1]
    [opt_value, opt_id] = min(J_temp,[],3);
    J_star = opt_value;
    if all(isnan(J_star(:)))
        disp('Table is too contaninated; no feasible solution');
    end
    
    savename=['COST_t=',num2str(t),'_Table.mat'];
    save(savename,'J_star','opt_id');
end
toc;
cd ..


disp('===============================================================');
disp('Starging Step 3');

%% initial condition
soc_ini = 0.55
We_ini = 0;

target_folder = [cyc_name, '_table_folder'];
cd(target_folder);

sim_We = zeros(1, length(cyc_time));
sim_Te = zeros(1, length(cyc_time));
sim_Wmg1 = zeros(1, length(cyc_time));
sim_Tmg1 = zeros(1, length(cyc_time));
sim_Wmg2 = zeros(1, length(cyc_time));
sim_Tmg2 = zeros(1, length(cyc_time));
sim_soc = zeros(1, length(cyc_time));
sim_fuel = zeros(1, length(cyc_time));

soc_c = soc_ini;
We_c = We_ini;
tic;
for t = 1:1:length(cyc_time)
    Wo_c = Wo(t); % current output speed (rad/s)
    To_c = To(t);
    dWo_c = dWo(t);
    
    if PowerD(t) <= 0 % vehicle is stopped or braking
        Tmg1_c = 0;
        Te_c = -6; % engine drag, which slows down the engine speed
    else % PowerD(t) > 0, driving
        load(['COST_t=',num2str(t),'_Table.mat']);
        id_lookup = interp2(x2_SOC_grid,x1_We_grid, opt_id, soc_c,We_c,'nearest');
        Tmg1_c = u1_Tmg1_grid(id_lookup);
        Te_c = interp1(We_optbsfc, Te_optbsfc, We_c);
    end
    Tmg2_c = (((dWo_c - T2W(2,1)*Te_c - T2W(2,3).*Tmg1_c)./T2W(2,2) + To_c))/MR;
    
    %% speed
    Wmg1_c = (We_c*D(1) + Wo_c*D(2))/(-D(3));
    Wmg2_c = Wo_c*MR;
    
    Tmg2_max = interp1(m_map_spd,m_max_trq,abs(Wmg2_c)); % [1]x[1]
    Tmg2_c(Tmg2_c>Tmg2_max) = Tmg2_max;
    Tmg2_c(Tmg2_c<-Tmg2_max) = -Tmg2_max;
    
    %% update x1 (engine speed)
    We_n = (T2W(1,1)*Te_c + T2W(1,2).*(Tmg2_c*MR - To_c) + T2W(1,3).*Tmg1_c)*dt + We_c;
    We_n = max(We_n, 0);
    fuel_c = dt*0.5*(interp2(eng_consum_trq,eng_consum_spd,eng_consum_fuel,Te_c,We_c,'linear')+...
                     interp2(eng_consum_trq,eng_consum_spd,eng_consum_fuel,Te_c,We_n,'linear')');
    
    %% update x2 (battery SOC)
    eff_mg1 = interp2(mg1_map_trq, mg1_map_spd, mg1_eff_map, Tmg1_c, abs(Wmg1_c));
    eff_mg2 = interp2(mg2_map_trq, mg2_map_spd, mg2_eff_map, Tmg2_c, abs(Wmg2_c));
    eff_mg1(isnan(eff_mg1))=0.2; % to avoid nan
    eff_mg2(isnan(eff_mg2))=0.2; % to avoid nan
    
    Pmg1_c = (Wmg1_c*Tmg1_c)*(eff_mg1^-sign(Wmg1_c*Tmg1_c));
    Pmg2_c = (Wmg2_c*Tmg2_c)*(eff_mg2^-sign(Wmg2_c*Tmg2_c));
    
    Pbatt = Pmg1_c + Pmg2_c; % +: discharge; -: charge
    if Pbatt>=0
        Pbatt_max = interp1(ess_soc, ess_max_pwr_dis, soc_c);
        Pbatt(Pbatt>Pbatt_max) = Pbatt_max;
        rint_c = interp1(ess_soc,ess_r_dis,soc_c);
    else
        Pbatt_min = interp1(ess_soc, ess_max_pwr_chg, soc_c);
        Pbatt(Pbatt<-Pbatt_min) = -Pbatt_min;
        rint_c = interp1(ess_soc,ess_r_chg,soc_c);
    end
    voc_c = interp1(ess_soc,ess_voc,soc_c);
    soc_n = (-(voc_c-(voc_c.^2-4*Pbatt.*rint_c).^(1/2))./(2*rint_c)/(ess_cap_ah*3600))*dt + soc_c;
    soc_n(soc_n>soc_max) = soc_max;
    
    %% ==============================
    sim_We(t) = We_c;
    sim_Te(t) = Te_c;
    sim_Wmg1(t) = Wmg1_c;
    sim_Tmg1(t) = Tmg1_c;
    sim_Tmg2(t) = Tmg2_c;
    sim_Wmg2(t) = Wmg2_c;
    sim_fuel(t) = fuel_c;
    sim_soc(t) = soc_c;
    
    We_c = We_n;
    soc_c = soc_n;
end
toc;
sim_Te(sim_Te<0) = 0;
dSOC = sim_soc(end) - sim_soc(1)

total_fuel_gram = sum(sim_fuel);
total_distance_mile = sum(cyc_mph)/3600;
mpg = total_distance_mile/(total_fuel_gram/1000/gasoline_density*liter2gallon)

cd ..


%% ================================================================== %%
figure(1); clf;
set(gcf, 'units', 'inch', 'pos', [7.7604    1.1146    5.8333    8.2813]);
subplot(4,1,1);
[haxes,hline1,hline2] = plotyy(cyc_time, cyc_mph(cyc_time), cyc_time,sim_soc);
set(hline1,'marker','x', 'markersize', 3, 'markerf', 'b');
set(hline2,'marker','x', 'markersize', 3, 'markerf', [0 0.5 0]);
set(haxes, 'fontsize', 8);
set(haxes(1), 'ylim', [0 80], 'ytick', 0:20:80, 'box', 'off', 'yminortick', 'on');
set(haxes(2), 'ylim', [0.4 0.7], 'ytick', 0.4:0.1:0.7, 'box', 'off', 'yminortick', 'on');
ylabel('Speed (mph)', 'parent', haxes(1));
ylabel('SOC (-)', 'parent', haxes(2));
set(haxes(2), 'XAxisLocation','top', 'xtick', [], 'xcolor', 'w');
set(gca, 'xgrid', 'on');
title(['\alpha = ', num2str(alpha), '; SOC_{ini} = ', num2str(sim_soc(1))]);

%% ==============================
subplot(4,1,2);
plot(cyc_time, sim_fuel, 'kx-', 'markersize', 3);
ylim([0 5]);
set(gca, 'pos', [0.1282    0.5935    0.7750    0.1495]);
set(gca, 'fontsize', 8);
set(gca, 'xgrid', 'on');
ylabel('Fuel Rate(g/s)');

%% ==============================
subplot(4,1,3);
plot(cyc_time, sim_We*rads2rpm, 'kx-', ...
     cyc_time, sim_Wmg1*rads2rpm, 'rx-', ...
     cyc_time, sim_Wmg2*rads2rpm, 'mx-', 'markersize', 3);
set(gca, 'pos', [0.1282    0.4171    0.7750    0.1495]);
ylim([-5000 15000]);
set(gca, 'fontsize', 8);
set(gca, 'xgrid', 'on');
ylabel('Speed (rpm)');
legend('Engine', 'MG1', 'MG2');

%% ==============================
subplot(4,1,4);
plot(cyc_time, sim_Te, 'kx-', ...
     cyc_time, sim_Tmg1, 'rx-', ...
     cyc_time, sim_Tmg2, 'mx-', 'markersize', 3);
set(gca, 'pos', [0.1282    0.2408    0.7750    0.1495]);
ylim([-100 150]);
set(gca, 'fontsize', 8);
set(gca, 'xgrid', 'on');
ylabel('Torque (Nm)');

% %% ==============================
% figure(2); clf;
% [C,h] = contour(eng_consum_spd*rads2rpm, eng_consum_trq, eng_bsfc', [216 200:5:250, 275, 300:100:800]);
% clabel(C,h, 'fontsize', 7);
% 
% hold on;
% plot(eng_map_spd*rads2rpm, eng_max_trq, 'k', 'linewidth', 2);
% plot(We_optbsfc*rads2rpm, Te_optbsfc, 'r--', 'linewidth', 2);
% plot(sim_We*rads2rpm, sim_Te, 'xb', 'markersize', 9, 'linewidth', 2);

