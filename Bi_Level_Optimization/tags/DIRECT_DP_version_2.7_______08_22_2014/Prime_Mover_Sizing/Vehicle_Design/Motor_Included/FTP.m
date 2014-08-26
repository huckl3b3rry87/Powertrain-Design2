% This program initially sizes the prime movers.

clc; clear all; close all;
global m rwh ni g grade Cd Af rho Frr eng_max_trq eng_consum_spd W_eng_min W_eng_max FD G
global m_map_spd m_max_trq Wm_max Wm_min nt

% Select a vehicle and components!
cd('Components');
TEST_Run_Sizing;  % If you change this the final solution changes..
Engine_41_kW;
% Engine_102_kW;
Battery_ADVISOR;
Motor_75_kW;
Vehicle_Parameters_4_HI_AV;
% Vehicle_Parameters_4_HI;
cd ..

% Desired Grade and Top Speed Performance - (Engine) - For Low Grades
% Test 1
V_eng = 50*mph_mps;
alpha_eng = 5*pi/180;
Pbase1 = Engine_Base_Power(V_eng,alpha_eng);

% Test 2
V_eng = 75*mph_mps;
alpha_eng = 0*pi/180;
Pbase2 = Engine_Base_Power(V_eng, alpha_eng);

% Initial Enigne Power
Pbase_kW =  max([Pbase1;Pbase2]);

% Calcualte the fc_trq_scale
fc_trq_scale = Pbase_kW/fc_max_pwr;  % For now neglect the effect that changing the mass of the engine has on this

% Desired Acceleration Performance
V_desired = 75.3;  % MPH
t_desired =  1; % Sec

% 72.1 to 74.9 in one second

% Get things started
Final_Velocity = 0;
n = 1;
mc_trq_scale =  1;

while(abs(Final_Velocity - V_desired) > 0.1)
    close all;
    if n > 1
        % Modify the motor torque scale based off actual accelaration performance
        mc_trq_scale = mc_trq_scale + (V_desired - Final_Velocity)/V_desired
    end
    % Enter proper vehicle parameters:
    cd('Components');
    Engine_41_kW;
    % Engine_102_kW;
    Motor_75_kW;                   % mass of battery will also change.. can do later..
    Vehicle_Parameters_4_HI_AV; % mass of motor will change in here..
    cd ..
    
    % Initialize Stuff
    V_sim = [];
    time_sim = [];
    
    % Gearbox and final drive ratios:
    ng= [4.484 2.872 1.842 1.414 1.000 0.742];
    
    x0=[72.1*mph_mps 0]; t0=0; tf = 0.5; % Initial conditions
    wem=eng_consum_spd(end-1); % Engine speed at which gear is shifted
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
            we(we < W_eng_min) = W_eng_min;
            we(we > W_eng_max) = W_eng_max;
            T_eng_max = interp1(eng_consum_spd,eng_max_trq,we);
            
            Wm_c = v/rwh*FD*G;
            Wm_c(Wm_c < Wm_min) = Wm_min;
            Wm_c(Wm_max < Wm_c) = Wm_max;
            Tm_max = interp1(m_map_spd,m_max_trq,Wm_c);
            
            Fti = T_eng_max*ni/rwh + Tm_max*FD*G/rwh;
%             Fti = Tm_max*FD*G/rwh;                         % Motor Only
            Frl = m*g*sin(grade) + Frr*m*g*cos(grade) + 0.5*rho*Cd*v.^2*Af;
            acc=(Fti-Frl)/m;% Differential equation for speed
            
            check = check + 1;

            if mean(acc) < 0.05
                break;
            end
            
            if we < wem
                tf = tf + 0.5;
            end
        
        end
        
        if check > 1  % Otherwise, just shift
            
            we_sim = 30*ni*v/rwh/pi; % RPM
            figure(3)
            plot(t, we_sim);
            hold on
            
            figure(4);
            plot(t,T_eng_max);
            hold on
            
            figure(5);
            plot(t,Tm_max);
            hold on
            
            figure(6);
            plot(t,Wm_c*30/pi);
            hold on
            
            figure(7);
            plot(t,Fti);
            hold on
            plot(t,Frl);
            hold on
            
            % At this point results for gear ‘ni’ are ready to plot
            figure(1);
            subplot(2,1,1), plot(t,v/mph_mps), hold on
            subplot(2,1,2), plot(t,s), hold on
            
            % Save Velocity and time
            V_sim = [V_sim; v(2:end)/mph_mps]; % MPH
            time_sim = [time_sim; t(2:end)];
            
            % Now gearshift is needed. Set the initial conditions for the next gear
            t0 = max(t);
            tf = tf + 0.5;
            x0=[v(end) s(end)];
            
            figure(2)
            plot(t, acc')
            hold on
            clear acc
        end
    end
    
    if check > 1
        
        figure(1)
        subplot(2,1,1)
        xlabel('Time (s)')
        ylabel('Velocity MPH')
        grid
        
        subplot(2,1,2)
        xlabel('Time (s)')
        ylabel('Distance (m)')
        grid
        
        figure(2)
        xlabel('Time (s)')
        ylabel('Acceleration (m/s^2)')
        grid
        
        figure(3)
        xlabel('Time (s)')
        ylabel('Engine Speed (RPM)')
        grid
        
        figure(4)
        xlabel('Time (s)')
        ylabel('Engine Torque (Nm)')
        grid
        
        figure(5)
        xlabel('Time (s)')
        ylabel('Motor Torque (Nm)')
        grid
        
        figure(6)
        xlabel('Time (s)')
        ylabel('Motor Speed (RPM)')
        grid
        
        figure(7);
        xlabel('Time (s)')
        ylabel('Force (N)')
        legend('Tractive Force','Road Load')
        grid
        
        figure(8)
        plot(time_sim,V_sim)
        xlabel('Time (s)')
        ylabel('Velocity MPH'), grid
        
        dt = 1;   % Careful with dt!
        [A,I] = min(abs(time_sim - dt)); 
        Final_Velocity = V_sim(I);
        
        if isempty(Final_Velocity)
            Final_Velocity = 0;
        end
        n = 2;
    end
end

% Add in ancilary power for both the engine and the motor
Paux_kW = Paux/1000;
mc_trq_scale = mc_trq_scale + Paux_kW/mc_max_pwr_initial_kW;
fc_trq_scale = fc_trq_scale + Paux_kW/fc_max_pwr_initial_kW;

% Recalcualte mc_max_pwr_kW & fc_mac_pwr_kW
cd('Components');
Engine_41_kW;
% Engine_102_kW;
Motor_75_kW;
cd ..

% At the lowest SOC for conservancy
Module_Number = 4*mc_max_pwr_kW*1000*Rint_size/(Voc_size^2)

mc_max_pwr_kW
mc_trq_scale
fc_max_pwr
fc_trq_scale