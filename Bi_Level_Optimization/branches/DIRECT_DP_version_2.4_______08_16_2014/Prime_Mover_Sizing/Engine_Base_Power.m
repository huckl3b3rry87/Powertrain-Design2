function [ Pbase ] = Engine_Base_Power(V_eng, alpha_eng)

global m g Frr rho Cd Af nt

 Pbase = (m*g*Frr*cos(alpha_eng) + 0.5*rho*Cd*Af*V_eng^2 + m*g*sin(alpha_eng))*V_eng/(1000*nt);
 

end

