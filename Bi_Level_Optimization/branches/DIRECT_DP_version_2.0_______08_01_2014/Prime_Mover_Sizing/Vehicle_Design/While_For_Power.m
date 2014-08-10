% Example 3.5.4
% Power estimation using numerical integration technique based on Figure 3.39
% Vehicle information:
m=1000; % Vehicle mass (kg)
c=0.4; % Aerodynamics overall coefficient
F0=200; % Rolling resistance force
theta=0; % Angle of slope (rad)
td=10; % Desired time (s)
vd=100; % Desired speed (km/h)

v0=eps; % Initial condition
t0=0; tf=10; % define differentiation interval t0-tf
% Initial value for Power:
Ps=(F0+c*vd^2)*vd+1;
% Iteration loop:
maxv=0;
% Continue iteration until the speed is very close to desired value ‘vd’
while abs(maxv-vd) > 0.001 % The answer is acceptable within a small tolerance
% Now invoke ode45:
[t,v]=ode45(@const_pow, [t0 tf], v0); % Calls ‘const_pow’ function
maxv=max(v);
% Check for velocity:
if maxv < vd % then increase power to increase velocity
P=P*vd/maxv;
else % decrease power (with a different rate) to decrease velocity
P=0.9*P*maxv/vd;
end
end
% At this point results are acceptable, so print them:
maxv
P_kW = P/1000