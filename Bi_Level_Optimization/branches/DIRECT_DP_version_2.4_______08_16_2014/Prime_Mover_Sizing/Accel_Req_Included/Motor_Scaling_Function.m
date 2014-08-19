function [ mc_trq_scale, Module_Number,Sim_Variables ] = Motor_Scaling_Function(V_0,V_desired, fc_trq_scale, dt_2)
global m rwh ni g grade Cd Af rho Frr eng_max_trq eng_consum_spd W_eng_min W_eng_max FD G
global m_map_spd m_max_trq Wm_max Wm_min nt Voc_size Rint_size

% Get things started
Final_Velocity = 1000;
n = 1;
mc_trq_scale =  1;

delta = abs(V_0 - V_desired)*0.1; % Within 5 percent
% delta = 0.1;
while(abs(Final_Velocity - V_desired) > delta)
    close all;
    if n > 1
        % Modify the motor torque scale based off actual accelaration performance
        mc_trq_scale = mc_trq_scale + (V_desired - Final_Velocity)/V_desired
    end
    
    % Enter proper vehicle parameters:
    cd('Components');
    
    % Eninge
    Engine_41_kW;
    % Engine_102_kW;
    
    % Motor
    Motor_75_kW;                 
   
    % Update the Battery Module Number based on The Current Motor Size
    module_number = 4*mc_max_pwr_kW*1000*Rint_size/(Voc_size^2);
    Battery_ADVISOR;
    
    % Calcualtes New Vehicle Mass
    Vehicle_Parameters_4_HI_AV; % mass of motor will change in here..
    cd ..
    
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
    wem=eng_consum_spd(end-2); % Engine speed at which gear is shifted
    we=0; % A value set for the ‘while’ statement below
    
    for i=1: length(ng) % Loop for gears
        ni=ng(i)*FD; % Overall gear ratio at each gear
        check = 0;
        we = 0;
        while we < wem % Repeat statements until we=wem
            [t,x]=ode45(@Fixed_thrt, [t0 tf], x0); % Calls Fixed_thrt function
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
end
Sim_Variables = [x_sim,V_sim,acc_sim,Teng_sim,Tm_sim,We_sim,Wm_sim,time_sim];

% Add in ancilary power for the motor
Paux_kW = Paux/1000;
mc_trq_scale = mc_trq_scale + Paux_kW/mc_max_pwr_initial_kW;

% Recalcualte mc_max_pwr_kW 
cd('Components');
Motor_75_kW;
cd ..

% At the lowest SOC for conservancy
Module_Number = 4*mc_max_pwr_kW*1000*Rint_size/(Voc_size^2);

end

