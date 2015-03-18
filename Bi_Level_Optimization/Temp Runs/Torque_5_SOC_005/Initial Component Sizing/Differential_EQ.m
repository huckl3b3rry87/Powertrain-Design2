function f = Differential_EQ(~, x, param, vinf, dvar, ni)

v = x(1); % instantaneous vehicle speed

we = ni*v/vinf.rwh; % rad/sec          % Does not depend on gear if you are only doing the motor!!
% we_sim = ni*v/vinf.rwh;
we(we < vinf.W_eng_min) = vinf.W_eng_min;
we(we > vinf.W_eng_max) = vinf.W_eng_max;
T_eng_max = interp1(vinf.eng_consum_spd_old,vinf.eng_max_trq,we);

wm = v/vinf.rwh*dvar.FD*dvar.G;
wm(wm < vinf.Wm_min) = vinf.Wm_min;
wm(vinf.Wm_max < wm) = vinf.Wm_max;
Tm_max = interp1(vinf.m_map_spd,vinf.m_max_trq,wm);

Fti = T_eng_max*ni/vinf.rwh + Tm_max*dvar.FD*dvar.G/vinf.rwh;
%             Fti = vinf.Tm_max*dvar.FD*dvar.G/vinf.rwh;                         % Motor Only
Frl = vinf.m*param.g*sin(param.grade) + vinf.Frr*vinf.m*param.g*cos(param.grade) + 0.5*param.rho*vinf.Cd*v.^2*vinf.Af;

acc = (Fti-Frl)/vinf.m;             % Differential equation for speed
dx_dt = v;               % Differential equation for distance  = dx_dt

f=[acc  
   dx_dt];