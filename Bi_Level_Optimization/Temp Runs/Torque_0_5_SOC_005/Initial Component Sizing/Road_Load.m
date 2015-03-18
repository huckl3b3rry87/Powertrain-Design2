function [ Road_Load ] = Road_Load(V, alpha_test, param, vinf )

Road_Load = (vinf.m*param.g*vinf.Frr*cos(alpha_test)*ones(size(V)) + 0.5*param.rho*vinf.Cd*vinf.Af*V.^2 + vinf.m*param.g*sin(alpha_test*ones(size(V))));


end

