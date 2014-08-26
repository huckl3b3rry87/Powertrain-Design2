function [ PASS, Sim_Variables ] = Acceleration_Test(V_0,V_f, Acc_Final, dt_2, param, vinf, dvar, TYPE)

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
time_sim =[];
in = 0;

% Gearbox and final drive ratios:
ng = vinf.gear;

x0=[V_0*param.mph_mps 0]; t0=0; tf = 0.25; % Initial conditions
wem = vinf.eng_consum_spd(end); % Engine speed at which gear is shifted

for g=1:1:length(ng) % Loop for gears
  
    ni=ng(g)*dvar.FD; % Overall gear ratio at each gear
    we = 0;
    while max(we) < wem % Repeat statements until we=wem
        [t,x]=ode23(@(t,x) Differential_EQ(t, x, param, vinf, dvar, ni), [t0 tf], x0); % Don't need to manipulate data again
        v=x(:,1);  % Speed
        s=x(:,2);  % Distance
        
        we = ni*v/vinf.rwh; % rad/sec          % Does not depend on gear if you are only doing the motor!!

        we(we < vinf.W_eng_min) = vinf.W_eng_min;
        
        we(we > vinf.W_eng_max) = vinf.W_eng_max;
        
        T_eng_max = interp1(vinf.eng_consum_spd,vinf.eng_max_trq,we);
        T_eng_max(we < 800*param.rpm2rads) = 0;
        
        wm = v/vinf.rwh*dvar.FD*dvar.G;
        wm(wm < vinf.Wm_min) = vinf.Wm_min;
        wm(vinf.Wm_max < wm) = vinf.Wm_max;
        Tm_max = interp1(vinf.m_map_spd,vinf.m_max_trq,wm);
        
        Fti = T_eng_max*ni/vinf.rwh + Tm_max*dvar.FD*dvar.G/vinf.rwh;
        %             Fti = vinf.Tm_max*dvar.FD*dvar.G/vinf.rwh;                         % Motor Only
        Frl = vinf.m*param.g*sin(param.grade) + vinf.Frr*vinf.m*param.g*cos(param.grade) + 0.5*param.rho*vinf.Cd*v.^2*vinf.Af;
        acc=(Fti-Frl)/vinf.m;             % Differential equation for speed
        
        if mean(acc) < 0.05  % Should hit the highest accelerations in the lowest gears
            break;  
        end
        
        if (max(we) < wem)
            tf = tf + 0.5;
        end
        
        if tf > 150
            break;    
        end
    end
    
    if any(we < wem) % Otherwise, the vehicle was traveling too fast, just shift
        % Save Variables
        V_sim = [V_sim; v(2:end)/param.mph_mps]; % MPH
        Fti_sim = [Fti_sim; Fti(2:end)];
        Frl_sim = [ Frl_sim;  Frl(2:end)];
        We_sim = [We_sim; we(2:end)*30/pi];   % RPM
        Wm_sim = [ Wm_sim; wm(2:end)*30/pi];   % RPM
        x_sim = [ x_sim; s(2:end)];
        acc_sim = [ acc_sim; acc(2:end)];
        Teng_sim = [Teng_sim; T_eng_max(2:end)];
        Tm_sim = [Tm_sim; Tm_max(2:end)];
        time_sim = [time_sim; t(2:end)];
        
        % Now gearshift is needed. Set the initial conditions for the next gear
        t0 = max(t);
        tf = tf + 1;
        x0=[v(end) s(end)];
        in = 1;
    end
end

if in == 1
    if TYPE == 1  % Final Velocity Req.
        [A,I] = min(abs(time_sim - dt_2));
        Final_Velocity = V_sim(I);
        
        if isempty(Final_Velocity)
            Final_Velocity = 0;
        end
        
        if Final_Velocity < V_f
            PASS = 0;            % Failed
        else
            PASS = 1;
        end
        
    else   % Acceleration req.
        Acc_Test = max(acc_sim);
        
        if isempty(Acc_Test)
            Acc_Test = 0;
        end
        
        if Acc_Test < Acc_Final
            PASS = 0;            % Failed
        else
            PASS = 1;
        end
        
    end
else
    PASS = 0;                   % Failed - Did not simulate for more than 0.5 s
end
Sim_Variables = [x_sim,V_sim,acc_sim,Teng_sim,Tm_sim,We_sim,Wm_sim,time_sim];
end





