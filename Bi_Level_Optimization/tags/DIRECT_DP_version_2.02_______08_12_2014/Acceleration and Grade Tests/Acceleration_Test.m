function [ PASS, Sim_Variables ] = Acceleration_Test(V_0,V_f, dt_2)

% Initialize Stuff
V_sim = [];
Fti_sim = [];
Frl_sim = [];
We_sim = [];
Wm_sim = [];
x_sim = [];
acc_sim = [];
Teng_sim = [];
Tm_sim = [];
time_sim = [];

% Gearbox and final drive ratios:
ng= [4.484 2.872 1.842 1.414 1.000 0.742];

x0=[V_0*mph_mps 0]; t0=0; tf = 0.5; % Initial conditions
wem = eng_consum_spd(end-2); % Engine speed at which gear is shifted
we = 0; % A value set for the ‘while’ statement below

for i=1: length(ng) % Loop for gears
    ni=ng(i)*FD; % Overall gear ratio at each gear
    check = 0;
    we = 0;
    while we < wem % Repeat statements until we=wem
        [t,x]=ode45(@Fixed_thrt, [t0 tf], x0);
        v=x(:,1);  % Speed
        s=x(:,2);  % Distance
        
        we = ni*v/rwh; % rad/sec          % Does not depend on gear if you are only doing the motor!!
        we_sim = ni*v/rwh;
        we(we < W_eng_min) = W_eng_min;
        we(we > W_eng_max) = W_eng_max;
        T_eng_max = interp1(eng_consum_spd,eng_max_trq,we);
        
        wm = v/rwh*FD*G;
        wm(wm < Wm_min) = Wm_min;
        wm(Wm_max < wm) = Wm_max;
        Tm_max = interp1(m_map_spd,m_max_trq,wm);
        
        Fti = T_eng_max*ni/rwh + Tm_max*FD*G/rwh;
        %             Fti = Tm_max*FD*G/rwh;                         % Motor Only
        Frl = m*g*sin(grade) + Frr*m*g*cos(grade) + 0.5*rho*Cd*v.^2*Af;
        acc=(Fti-Frl)/m;             % Differential equation for speed
        
        check = check + 1;
        
        if mean(acc) < 0.05
            break;
        end
        
        if (we < wem)
            tf = tf + 0.5;
        end
        
    end
    
    if check > 1  % Otherwise, just shift
        
        % Save Variables
        V_sim = [V_sim; v(2:end)/mph_mps]; % MPH
        Fti_sim = [Fti_sim; Fti(2:end)];
        Frl_sim = [ Frl_sim;  Frl(2:end)];
        We_sim = [We_sim; we_sim(2:end)*30/pi];   % RPM
        Wm_sim = [ Wm_sim; wm(2:end)*30/pi];   % RPM
        x_sim = [ x_sim; s(2:end)];
        acc_sim = [ acc_sim; acc(2:end)];
        Teng_sim = [Teng_sim; T_eng_max(2:end)];
        Tm_sim = [Tm_sim; Tm_max(2:end)];
        time_sim = [time_sim; t(2:end)];
        
        % Now gearshift is needed. Set the initial conditions for the next gear
        t0 = max(t);
        tf = tf + 0.5;
        x0=[v(end) s(end)];
        
        clear acc
    end
end

if check > 1
    
    [A,I] = min(abs(time_sim - dt_2));
    Final_Velocity = V_sim(I);
    
    if isempty(Final_Velocity)
        Final_Velocity = 0;
    end
    n = 2;
end

Sim_Variables = [x_sim,V_sim,acc_sim,Teng_sim,Tm_sim,We_sim,Wm_sim,time_sim];

if Final_Velocity < V_f
    PASS = 0;
else
    PASS = 1;
end
end





