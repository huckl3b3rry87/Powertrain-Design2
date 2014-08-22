function [ P_eng_mass ] = Engine_Mass_Power(V_eng, alpha_eng, param, vinf, delta_fc_mass)

P_eng_mass = (delta_fc_mass*param.g*vinf.Frr*cos(alpha_eng) + delta_fc_mass*param.g*sin(alpha_eng))*V_eng/(1000*vinf.nt);  % Only add or subrtact the things that change from mass!!

end

