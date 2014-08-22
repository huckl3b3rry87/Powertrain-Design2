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
time_sim = [];

% Gearbox and final drive ratios:
ng= [4.484 2.872 1.842 1.414 1.000 0.742];

x0=[V_0*param.mph_mps 0]; t0=0; tf = 0.25; % Initial conditions
wem = vinf.eng_consum_spd(end-2); % Engine speed at which gear is shifted


% global Fti Frl we_sim wm T_eng_max Tm_max acc 

for i=1: length(ng) % Loop for gears
    
    if i > 1
        if max(v)/param.mph_mps > V_f*1.3
            break;
        end
    end
    
    ni=ng(i)*dvar.FD; % Overall gear ratio at each gear
    check = 0;
    we = 0;
    while we < wem % Repeat statements until we=wem
        [t,x]=ode23(@(t,x) Differential_EQ(t, x, param, vinf, dvar, ni), [t0 tf], x0); % Don't need to manipulate data again
        v=x(:,1);  % Speed
        s=x(:,2);  % Distance
        
%         f = Differential_EQ(t, x, param, vinf, dvar, ni);
%         acc = f(2);
        
        we = ni*v/vinf.rwh; % rad/sec          % Does not depend on gear if you are only doing the motor!!
        we_sim = ni*v/vinf.rwh;
        we(we < vinf.W_eng_min) = vinf.W_eng_min;
        we(we > vinf.W_eng_max) = vinf.W_eng_max;
        T_eng_max = interp1(vinf.eng_consum_spd,vinf.eng_max_trq,we);
        
        wm = v/vinf.rwh*dvar.FD*dvar.G;
        wm(wm < vinf.Wm_min) = vinf.Wm_min;
        wm(vinf.Wm_max < wm) = vinf.Wm_max;
        Tm_max = interp1(vinf.m_map_spd,vinf.m_max_trq,wm);
        
        Fti = T_eng_max*ni/vinf.rwh + Tm_max*dvar.FD*dvar.G/vinf.rwh;
        %             Fti = vinf.Tm_max*dvar.FD*dvar.G/vinf.rwh;                         % Motor Only
        Frl = vinf.m*param.g*sin(param.grade) + vinf.Frr*vinf.m*param.g*cos(param.grade) + 0.5*param.rho*vinf.Cd*v.^2*vinf.Af;
        acc=(Fti-Frl)/vinf.m;             % Differential equation for speed
        
        check = check + 1;
        if mean(acc) < 0.05
            TOO_LONG = 0;
            break; 
        end
        
        if (we < wem)
            tf = tf + 1;
        end
        
        if tf > 100
            TOO_LONG = 1;
            break;    
        else
            TOO_LONG = 0;
        end
    end
    
    if check > 1  % Otherwise, just shift
        % Save Variables
        V_sim = [V_sim; v(2:end)/param.mph_mps]; % MPH
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

if check > 1 && TOO_LONG ~= 1
    
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





