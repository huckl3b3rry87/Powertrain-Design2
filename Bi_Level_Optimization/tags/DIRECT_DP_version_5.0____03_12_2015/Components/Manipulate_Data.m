% Motor
m_map_trq = mc_trq_scale*m_map_trq;
m_max_trq = mc_trq_scale*m_max_trq;
m_max_gen_trq = mc_trq_scale*m_max_gen_trq;
mc_mass = mc_trq_scale*mc_mass;

% Engine
eng_control_trq = fc_trq_scale*eng_control_trq;
eng_consum_trq = fc_trq_scale*eng_consum_trq;
eng_max_trq = fc_trq_scale*eng_max_trq;
eng_consum_fuel = fc_trq_scale*eng_consum_fuel;
Engine_Mass;

% Battery
ess_r_dis = module_number*ess_r_dis;
ess_r_chg = module_number*ess_r_chg;
ess_voc = module_number*ess_voc;
ess_mass = module_number*ess_module_mass;
ess_max_pwr_dis = (ess_voc.^2)./(4*ess_r_dis)*0.98;  
ess_max_pwr_chg = (ess_voc.^2)./(4*ess_r_chg)*0.98;

% Calcualte Total Mass
m = mc_mass + fc_mass + Base_Vehicle + Load + ess_mass;