% Function called in program main_ft.m
function f=Fixed_thrt(t,x)
global m rw ni g grade Cd Af rho Frr eng_max_trq eng_consum_spd W_eng_min W_eng_max
v=x(1); % instantaneous vehicle speed
We_c = ni*v*30/rw/pi; % rad/sec

if We_c < W_eng_min
    We_c = W_eng_min;
end
if We_c > W_eng_max
    We_c = W_eng_max;
end

Trq = interp1(eng_consum_spd,eng_max_trq,We_c);
Fti=Trq*ni/rw;
Frl = m*g*sin(grade) + Frr*m*g*cos(grade) + 0.5*rho*Cd*(v)^2*Af;

f1=(Fti-Frl)/m;% Differential equation for speed = dv_dt
f2=v; % Differential equation for distance  = dx_dt

f=[f1 
   f2];