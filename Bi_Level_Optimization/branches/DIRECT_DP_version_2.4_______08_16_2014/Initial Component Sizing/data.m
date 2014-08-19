% Put all the data into a structure

% Universal Parameters
param.g = 9.81;   % m/s^2
param.rho = 1.2;  % density of air (kg/m^3)
param.mph_mps = 1/2.237;
param.gasoline_density = 0.7197; % [kg/liter]
param.liter2gallon = 0.264172;
param.MIN_SOC = 0.4;
param.MAX_SOC = 0.8;
param.grade = 0;

% Vehicle information
vinf.base_mass = Base_Vehicle;
vinf.load = Load;
vinf.rwh = rwh;
vinf.Frr = Frr;
vinf.Cd = Cd;
vinf.Af = Af;
vinf.Paux = Paux;
vinf.nt = nt;

% Motor
vinf.m_map_spd = m_map_spd;
vinf.m_map_trq_orig = m_map_trq;
vinf.m_max_trq_orig = m_max_trq;
vinf.m_max_gen_trq_orig = m_max_gen_trq;
vinf.m_eff_map = m_eff_map;
vinf.Wm_min = Wm_min;
vinf.Wm_max = Wm_max;
vinf.mc_mass_orig = mc_mass;
vinf.mc_max_pwr_kW = mc_max_pwr_kW;

% Engine
vinf.eng_control_trq_orig= eng_control_trq;
vinf.W_eng_min = W_eng_min;
vinf.W_eng_max = W_eng_max;
vinf.eng_consum_spd = eng_consum_spd;
vinf.eng_consum_trq_orig = eng_consum_trq;
vinf.eng_max_trq_orig = eng_max_trq;
vinf.optimal_eng_spd = optimal_eng_spd;
vinf.Te_min = Te_min;
vinf.fc_max_pwr = fc_max_pwr;
vinf.eng_consum_fuel_orig = eng_consum_fuel;

% Battery
vinf.ess_r_dis_orig = ess_r_dis;
vinf.ess_r_chg_orig = ess_r_chg;
vinf.ess_voc_orig = ess_voc;
vinf.ess_cap_ah = ess_cap_ah;
vinf.ess_soc = ess_soc;
vinf.ess_max_pwr_dis = ess_max_pwr_dis;
vinf.ess_max_pwr_chg = ess_max_pwr_chg;
vinf.ess_module_mass_orig = ess_module_mass;

