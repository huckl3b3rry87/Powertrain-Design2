
% ADVISOR data file:  ESS_LI7_temp.m
% 
%module_number = 63; % Could round later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temperature range over which data is defined
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ess_tmp=[0 25 41];  % (C)


ess_r_dis=[0.0419 0.0288 0.0221 0.014 0.0145 0.0145 0.0162;
0.072 0.01515 0.00839 0.00493 0.00505 0.005524 0.005722;
0.0535 0.0133 0.0082 0.0059 0.0059 0.006 0.0063]*3*module_number; % (ohm)
% module's resistance to being charged, indexed by ess_soc and ess_tmp
ess_r_chg=[0.021 0.018 0.0177 0.0157 0.0138 0.0138 0.015;
0.0124 0.0068 0.005426 0.00442 0.00463 0.00583 0.00583;
0.0104 0.0079 0.0072 0.0064 0.0059 0.0058 0.006]*3*module_number; % (ohm)
% module's open-circuit (a.k.a. no-load) voltage, indexed by ess_soc and ess_tmp
ess_voc=[3.44 3.473 3.496 3.568 3.637 3.757 3.896;
3.124 3.349 3.433 3.518 3.616 3.752 3.898;
3.128 3.36 3.44 3.528 3.623 3.761 3.899]*3*module_number; % (V)

ess_cap_ah=[5.943 7.035 7.405]; % (A*h), max. capacity at C/3 rate, indexed by ess_tmp

%Assunme that the battery pack is 41 C
ess_r_dis = ess_r_dis(3,:);
ess_r_chg = ess_r_chg(3,:);
ess_voc = ess_voc(3,:);
ess_cap_ah = ess_cap_ah(3);

% For Battery sizing
Rint_size = [ess_r_dis(4)/module_number];
Voc_size =  [ess_voc(4)/module_number];

ess_soc=[0 10 20 40 60 80 100]/100; 

% BORROWED the resistance down scaling helps to match the 100kW net power spec of 2010 Prius (mentioned in Namwook Kim's journal)
ess_max_pwr_dis = (ess_voc.^2)./(4*ess_r_dis)*0.98; % 110~105kW in 40-70 percent SOC window
ess_max_pwr_chg = (ess_voc.^2)./(4*ess_r_chg)*0.98; % 105~100kW in 40-70 percent SOC window

% plot(ess_soc,ess_max_pwr_dis/1000)

ess_module_mass = 0.37824*3;  % (kg), mass of single module

ess_mass = module_number*ess_module_mass;