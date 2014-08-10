% Function called in program main_ft.m
function f=Fixed_thrt(t,x)
global m rwh ni g grade Cd Af rho Frr eng_max_trq eng_consum_spd W_eng_min W_eng_max FD G
global m_map_spd m_max_trq Wm_max Wm_min
v=x(1); % instantaneous vehicle speed

% ENGINE
We_c = ni*v/rwh; % rad/sec
if We_c < W_eng_min
    We_c = W_eng_min;
end
if We_c > W_eng_max
    We_c = W_eng_max;
end
T_eng_max = interp1(eng_consum_spd,eng_max_trq,We_c);

% MOTOR
Wm_c = v/rwh*FD*G;
if Wm_c > Wm_max
    Wm_c = Wm_max;
end
if Wm_c < Wm_min
    Wm_c = Wm_min;
end
Tm_max = interp1(m_map_spd,m_max_trq,Wm_c);

% Total
Fti = T_eng_max*ni/rwh + Tm_max*FD*G/rwh;

Frl = m*g*sin(grade) + Frr*m*g*cos(grade) + 0.5*rho*Cd*v^2*Af;

f1 = (Fti-Frl)/m;% Differential equation for speed = dv_dt
f2 = v; % Differential equation for distance  = dx_dt

f=[f1 
   f2];