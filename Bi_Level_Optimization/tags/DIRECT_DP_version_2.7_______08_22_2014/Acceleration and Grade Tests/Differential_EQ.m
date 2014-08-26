function f = Differential_EQ(~, x, param, vinf, dvar, ni)

v = x(1); % instantaneous vehicle speed

% ENGINE
We_c = ni*v/vinf.rwh; % rad/sec
if We_c < vinf.W_eng_min
    We_c = vinf.W_eng_min;
end
if We_c > vinf.W_eng_max
    We_c = vinf.W_eng_max;
end
T_eng_max = interp1(vinf.eng_consum_spd,vinf.eng_max_trq,We_c);

% MOTOR
Wm_c = v/vinf.rwh*dvar.FD*dvar.G;
if Wm_c > vinf.Wm_max
    Wm_c = vinf.Wm_max;
end
if Wm_c < vinf.Wm_min
    Wm_c = vinf.Wm_min;
end
Tm_max = interp1(vinf.m_map_spd,vinf.m_max_trq,Wm_c);

% Total
Fti = T_eng_max*ni/vinf.rwh + Tm_max*dvar.FD*dvar.G/vinf.rwh;
% Fti = Tm_max*FD*G/rwh;                         % Motor Only

Frl = vinf.m*param.g*sin(param.grade) + vinf.Frr*vinf.m*param.g*cos(param.grade) + 0.5*param.rho*vinf.Cd*v^2*vinf.Af;

f1 = (Fti-Frl)/vinf.m;% Differential equation for speed = dv_dt
f2 = v;               % Differential equation for distance  = dx_dt

f=[f1 
   f2];