function [ Pbase ] = Engine_Base_Power(V_eng, alpha_eng, param, vinf)

Pbase = (vinf.m*param.g*vinf.Frr*cos(alpha_eng) + 0.5*param.rho*vinf.Cd*vinf.Af*V_eng^2 + vinf.m*param.g*sin(alpha_eng))*V_eng/(vinf.nt);

end

