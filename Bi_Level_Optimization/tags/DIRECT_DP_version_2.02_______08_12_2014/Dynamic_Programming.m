% Dynamic Programming - written by Huckleberry Febbo - 07/20/2014 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%----------------------------DEFINE THE RUN TYPE--------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

RUN_TYPE = 0;  % RUN_TYPE = 1 - for DIRECT     &    RUN_TYPE = 0 - for DP only

if RUN_TYPE == 0
    clear all
    close all
    clc
    RUN_TYPE = 0;
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%----------------------------Load All Data--------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

cd('Components');
if RUN_TYPE == 0
    %TEST_Run_4_HI;
    TEST_Run_4_HI_AV;
end
%                              ~~ Engine ~~
% Engine_2rz_0410; 
% Engine_102_kW;
Engine_41_kW;

%                              ~~ Motor ~~
% Motor_int;
Motor_75_kW;

%                             ~~ Battery ~~
% Battery_int;  % No variation with the number of modules in this battery!!
Battery_ADVISOR;

%                              ~~ Vehicle ~~

Vehicle_Parameters_4_HI_AV;
% Vehicle_Parameters_4_HI;
cd ..  

cd('Drive_Cycle')

%load CYC_HWFET; cyc_name = 'HWFET';
%load CYC_UDDS; cyc_name = 'UDDS';
load CYC_US06; cyc_name = 'US06';
%load SHORT_CYC_HWFET; cyc_name = 'SHORT_CYC_HWFET';
%load RAMP; cyc_name = 'RAMP';
%load RAMP_slow; cyc_name = 'RAMP_slow';
% load CYC_LA92; cyc_name = 'LA92';
%load CONST_65; cyc_name = 'CONST_65';
%load CONST_45; cyc_name = 'CONST_45';
%load CYC_COMMUTER; cyc_name = 'COMMUTER';

% AV Cycles
% load CYC_US06_AV; cyc_name = 'US06_AV';

Manipulate_Drive_Cycle;
cd ..
tables = [cyc_name, ' TABLES'];
mkdir(tables)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%---------------------Simulate all Possible Dynamics----------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%Define State Grids
x1_grid = [0.4:0.005:0.8]';       % SOC 
x1_length = length(x1_grid);

x2_grid = [0 1];   % Engine off (0)   &   Engine on (1)
x2_length = length(x2_grid);

x3_grid = [4.484 2.872 1.842 1.414 1.000 0.742];  % [1st 2nd...]             % Gear Level
x3_length = length(x3_grid);

% Define Control Grids
u1_grid = eng_control_trq';          % Engine Torque [N-m]
u1_length = length(u1_grid);  

u2_grid = [-1 0 1];             % Shift down, nothing, up
u2_length = length(u2_grid);

u3_length = 2;

%% Define Some soft cnstraints stuff
SOC_penalty = linspace(0.1,10,20);
NEAR_SOC_min = MIN_SOC + fliplr(linspace(0.001,0.02,20));
NEAR_SOC_max = MAX_SOC - fliplr(linspace(0.001,0.02,20));

% Make a New Folder 
rmdir(tables,'s')                 % Delete any left over info
tables = [cyc_name, ' TABLES'];
mkdir(tables)
cd(tables);

if RUN_TYPE == 0
    tic
end

for t = 1:time_cyc
    Tm_max = single(zeros(x2_length,x3_length,u1_length,u2_length,u3_length));               % [x2]x[u1]x[u2]
    Tm_save = single(zeros(x2_length,x3_length,u1_length,u2_length,u3_length));              % [x2]x[u1]x[u2]
    Wm_save = single(zeros(x2_length,x3_length,u1_length,u2_length,u3_length));              % [x2]x[u1]x[u2]
    We_save = single(zeros(x2_length,x3_length,u1_length,u2_length,u3_length));              % [x2]x[u1]x[u2]
    table_x1 = single(zeros(x1_length,x2_length,x3_length,u1_length,u2_length,u3_length));   % [x2]x[u1]x[u2]
    inst_fuel = single(zeros(x2_length,x3_length,u1_length,u2_length,u3_length));            % [x2]x[u1]x[u2]
    infeasible_Te = single(zeros(x2_length,x3_length,u1_length,u2_length,u3_length));        % [x2]x[u1]x[u2]
    infeasible_Pbatt = single(zeros(x1_length,x2_length,x3_length,u1_length,u2_length,u3_length));
    
    for x3 = 1:x3_length           % Go through all of the gears
        x3_c = x3_grid(x3);
        
        for u2 = 1:u2_length       % Shift down, don't shift, and shift up
            u2_c = u2_grid(u2);
            
            if (x3_c == x3_grid(1) && u2_c == u2_grid(1)) || (x3_c == x3_grid(x3_length) && u2_c == u2_grid(u2_length))
                u2_c = 0;          % Cannot shift
            end
            
            if u2 == 1 || u2 == 3  % Shift penalty
                Shift_Penalty = 0.25;
            else
                Shift_Penalty = 0;
            end
            
            % Update x3 and Define New Gear ID
            New_Gear_Index = x3 + u2_c;
            x3_n = x3_grid(New_Gear_Index);
            
            for x2 = 1:x2_length             % Current Engine State
                ENG_state_c = x2_grid(x2);
                
                for u3 = 1:u3_length         % Engine Control
                    
                    if  ENG_state_c == 0 && u3 == 2   % Engine will be turned on - Done later in the code..this is the same thing though.
                        Eng_Penalty = fc_trq_scale*20;
                    else
                        Eng_Penalty = 0;
                    end
                    
                    if ENG_state_c == 0 && u3 == 1 || ENG_state_c == 1 && u3 == 1    % Next Engine State is off
                        ENG_state_n = 0;
                        We_c = 0;                            % [rad/sec]
                        Wm_c = Ww(t)*FD*G;              % [rad/sec]
                        Te_c = zeros(u1_length,1);           % Turn the engine off
                    else                                                             % Engine is on
                        ENG_state_n = 1;
                        We_c = Ww(t)*FD*x3_n;                % [rad/sec]
                        Taux = Paux/We_c;
                        Te_c =  u1_grid - Taux;                      % Engine Control, [u1]x[1]- What if this becomes negative!!
                        Wm_c = Ww(t)*FD*G;              % [rad/sec]
                    end
                    
                    if Pd(t) < 0 && ENG_state_n == 1 % Braking and Running Engine - 
                        We_c = optimal_eng_spd;
                        Taux = Paux/We_c;
                        Te_c = Taux*ones(size(Te_c));                    
                        Tm_c = Tw(t)/(FD*G)*ones(size(Te_c)); % Engine is not connected to motor
                    else
                        Tm_c = Tw(t)/(FD*G)*ones(size(Te_c)) - Te_c*x3_n/G;  % [u1]x[1]
                    end
                    
                    % Check Motor
                    Tm_max_current = interp1(m_map_spd,m_max_trq,Wm_c)*ones(size(u1_grid));
                    Tm_max(x2,x3,:,u2,u3) =  Tm_max_current;    % [x2]x[x3]x[u1]x[u2]x[u3]   - #Check to see if the motor map is symmetric
                    Tm_save(x2,x3,:,u2,u3) = Tm_c;                                                      % [x2]x[x3]x[u1]x[u2]x[u3]
                    Wm_save(x2,x3,:,u2,u3) = Wm_c*ones(size(u1_grid));                                  % [u1]x[1]  -  Check For Each Gear
                    
                    % Check Engine
                    We_save(x2,x3,:,u2,u3) = We_c*ones(size(u1_grid));                         % [u1]x[1]  - Check For Each Gear - Speed of the engine depends on the gear and engine state
                   
                    % Saturate Engine Speed - For Max Torque Lookup & Fuel Lookup
                    if We_c < W_eng_min
                        We_c = W_eng_min;
                    end
                    if We_c > W_eng_max
                        We_c = W_eng_max;
                    end
                    
                    Te_max =  interp1(eng_map_spd,eng_max_trq,We_c)*ones(size(Te_c));
                    
                    if u3 == 2  % The Engine is on
                        infeasible_Te(x2,x3,:,u2,u3) = (Te_max < Te_c);                            % [u1]x[1]
                    else
                        infeasible_Te(x2,x3,:,u2,u3) = zeros(size(Te_c));
                    end
                    
                    % Saturate Engine Torque - For Fuel Lookup
                    Te_c(Te_c > Te_max) = Te_max(Te_c > Te_max);
                    Te_c(Te_c < Te_min) = Te_min(Te_c < Te_min);    % saturate it, but do not penalize it!
                    
                    if ENG_state_n == 1
                        fuel = (interp2(eng_consum_trq,eng_consum_spd,eng_consum_fuel,Te_c,We_c,'linear')*dt)';                 
                    else
                        fuel = zeros(size(Te_c));
                    end
                    
                    inst_fuel(x2,x3,:,u2,u3) = fuel + Shift_Penalty*ones(size(fuel)) + Eng_Penalty*ones(size(fuel));    
                     
                    % Update x1
                    
                    % Saturate the motor for the efficiency lookup table
                    Tm_c(Tm_c > Tm_max_current) = Tm_max_current(Tm_c > Tm_max_current);
                    Tm_c(Tm_c < -Tm_max_current) = -Tm_max_current(Tm_c < -Tm_max_current);    
                    Wm_c(Wm_c > Wm_max) = Wm_max;
                    Wm_c(Wm_c < Wm_min) = Wm_min;
                    
                    eff_m = interp2(m_map_trq, m_map_spd, m_eff_map, Tm_c, abs(Wm_c))';      
                    eff_m(isnan(eff_m)) = 0.2;
                    
                    Pbat_charge = (Wm_c*Tm_c).*(eff_m);    % Tm_c < 0
                    Pbat_discharge = (Wm_c*Tm_c)./(eff_m);
                    
                    Pbat = Pbat_discharge;
                    Pbat(Tm_c < 0) = Pbat_charge(Tm_c < 0);
                    
                    Pbat = repmat(Pbat,[1, x1_length]);                                      
                    Pbat = permute(Pbat, [2 1]);
                    
                    if ENG_state_n == 0
                        Pbat = Pbat + Paux;    % Draw Aux Power From Battery
                    end
                    
                    % Discharge
                    Pbatt_max = repmat(interp1(ess_soc, ess_max_pwr_dis, x1_grid),[1,u1_length]);
                    rint_discharge = repmat(interp1(ess_soc,ess_r_dis,x1_grid),[1,u1_length]);
                    
                    % Charge
                    Pbatt_min = -repmat(interp1(ess_soc, ess_max_pwr_chg, x1_grid),[1,u1_length]);
                    rint_charge = repmat(interp1(ess_soc,ess_r_chg,x1_grid),[1,u1_length]);
                  
                    % Check Battery Infeasibility
                    infeasible_Pbatt(:,x2,x3,:,u2,u3) = (Pbatt_max < Pbat); % Do not peanlize it for opperating too low, we can brake
                     
                    % Saturate the Battery
                    Pbat(Pbatt_max < Pbat) = Pbatt_max(Pbatt_max < Pbat);
                    Pbat(Pbat < Pbatt_min) = Pbatt_min(Pbat < Pbatt_min);
                    
                    % Charge & Discharge Resistances
                    rint_c = rint_charge;
                    rint_c(Pbat > 0) = rint_discharge(Pbat > 0);
                    
                    Voc_c = repmat(interp1(ess_soc,ess_voc,x1_grid), [1, u1_length]);
                    SOC_c_matrix = repmat(x1_grid,[1, u1_length]);
                    SOC_n =  SOC_c_matrix -(Voc_c -(Voc_c.^2 -4*Pbat.*rint_c).^(1/2))./(2*rint_c*ess_cap_ah*3600)*dt;
                    
                    table_x1(:,x2,x3,:,u2,u3) = SOC_n;  
                                      
                end   % End of u3 ( Engine Control ) Loops
            end       % End of x2 ( Engine State ) Loops
        end           % End of u2 ( Gear Control )
    end               % End of x3 (Gear State) loops
    % Check Motor
%     infeasible_Tm = (Tm_save > Tm_max) | (Tm_save < -Tm_max);             % [x2]x[x3]x[u1]x[u2]x[u3] 
    infeasible_Tm = (Tm_save > Tm_max);            % Can Brake to make the rest up
    infeasible_Tm = repmat(infeasible_Tm,[1,1,1,1,1,x1_length]);
    infeasible_Tm = permute(infeasible_Tm,[6 1 2 3 4 5]);
    
    infeasible_Wm = (Wm_save > Wm_max) | (Wm_save < Wm_min);              % [x2]x[x3]x[u1]x[u2]x[u3]
    infeasible_Wm = repmat(infeasible_Wm,[1,1,1,1,1,x1_length]);
    infeasible_Wm = permute(infeasible_Wm,[6 1 2 3 4 5]);
    
    % Check Engine 
    infeasible_We = single(zeros(x2_length,x3_length,u1_length,u2_length,u3_length));
    infeasible_We(:,:,:,:,2) = (We_save(:,:,:,:,2) < 800*rpm2rads) | (We_save(:,:,:,:,2) > W_eng_max);
    infeasible_We = repmat(infeasible_We,[1,1,1,1,1,x1_length]);
    infeasible_We = permute(infeasible_We,[6 1 2 3 4 5]);
    
    infeasible_Te = repmat(infeasible_Te,[1,1,1,1,1,x1_length]);
    infeasible_Te = permute(infeasible_Te,[6 1 2 3 4 5]);
    
    % Check SOC
    infeasible_SOC = (table_x1 < MIN_SOC) | (table_x1 > MAX_SOC);        % [x2]x[x3]x[u1]x[u2]x[u3]
    table_x1(table_x1 > MAX_SOC) = MAX_SOC; 
    table_x1(table_x1 < MIN_SOC) = MIN_SOC;
    
    % Lower SOC penalties
    SOC_soft = SOC_penalty(1)*((NEAR_SOC_min(2) < table_x1)  & (table_x1 < NEAR_SOC_min(1))); % Check this later!!
    SOC_soft = SOC_soft + SOC_penalty(2)*(NEAR_SOC_min(3) < table_x1  & table_x1 < NEAR_SOC_min(2));
    SOC_soft = SOC_soft + SOC_penalty(3)*(NEAR_SOC_min(4) < table_x1  & table_x1 < NEAR_SOC_min(3));
    SOC_soft = SOC_soft + SOC_penalty(4)*(NEAR_SOC_min(5) < table_x1  & table_x1 < NEAR_SOC_min(4));
    SOC_soft = SOC_soft + SOC_penalty(5)*(NEAR_SOC_min(6) < table_x1  & table_x1 < NEAR_SOC_min(5));
    SOC_soft = SOC_soft + SOC_penalty(6)*(NEAR_SOC_min(7) < table_x1  & table_x1 < NEAR_SOC_min(6));
    SOC_soft = SOC_soft + SOC_penalty(7)*(NEAR_SOC_min(8) < table_x1  & table_x1 < NEAR_SOC_min(7));
    SOC_soft = SOC_soft + SOC_penalty(8)*(NEAR_SOC_min(9) < table_x1  & table_x1 < NEAR_SOC_min(8));
    SOC_soft = SOC_soft + SOC_penalty(9)*(NEAR_SOC_min(10)< table_x1  & table_x1 < NEAR_SOC_min(9));
    SOC_soft = SOC_soft + SOC_penalty(10)*(NEAR_SOC_min(11)< table_x1 & table_x1 < NEAR_SOC_min(10));
    SOC_soft = SOC_soft + SOC_penalty(11)*(NEAR_SOC_min(12)< table_x1 & table_x1 < NEAR_SOC_min(11));
    SOC_soft = SOC_soft + SOC_penalty(12)*(NEAR_SOC_min(13)< table_x1 & table_x1 < NEAR_SOC_min(12));
    SOC_soft = SOC_soft + SOC_penalty(13)*(NEAR_SOC_min(14)< table_x1 & table_x1 < NEAR_SOC_min(13));
    SOC_soft = SOC_soft + SOC_penalty(14)*(NEAR_SOC_min(15)< table_x1 & table_x1 < NEAR_SOC_min(14));
    SOC_soft = SOC_soft + SOC_penalty(15)*(NEAR_SOC_min(16)< table_x1 & table_x1 < NEAR_SOC_min(15));
    SOC_soft = SOC_soft + SOC_penalty(16)*(NEAR_SOC_min(17)< table_x1 & table_x1 < NEAR_SOC_min(16));
    SOC_soft = SOC_soft + SOC_penalty(17)*(NEAR_SOC_min(18)< table_x1 & table_x1 < NEAR_SOC_min(17));
    SOC_soft = SOC_soft + SOC_penalty(18)*(NEAR_SOC_min(19)< table_x1 & table_x1 < NEAR_SOC_min(18));
    SOC_soft = SOC_soft + SOC_penalty(19)*(NEAR_SOC_min(20)< table_x1 & table_x1 < NEAR_SOC_min(19));
    SOC_soft = SOC_soft + SOC_penalty(20)*((MIN_SOC < table_x1) & (table_x1 < NEAR_SOC_min(20)));
    
    % Upper SOC penalties
    SOC_soft = SOC_soft + SOC_penalty(1)*(NEAR_SOC_max(2)  > table_x1 & table_x1 > NEAR_SOC_max(1)); % Check this later!!
    SOC_soft = SOC_soft + SOC_penalty(2)*(NEAR_SOC_max(3)  > table_x1 & table_x1 > NEAR_SOC_max(2));
    SOC_soft = SOC_soft + SOC_penalty(3)*(NEAR_SOC_max(4)  > table_x1 & table_x1 > NEAR_SOC_max(3));
    SOC_soft = SOC_soft + SOC_penalty(4)*(NEAR_SOC_max(5)  > table_x1 & table_x1 > NEAR_SOC_max(4));
    SOC_soft = SOC_soft + SOC_penalty(5)*(NEAR_SOC_max(6)  > table_x1 & table_x1 > NEAR_SOC_max(5));
    SOC_soft = SOC_soft + SOC_penalty(6)*(NEAR_SOC_max(7)  > table_x1 & table_x1 > NEAR_SOC_max(6));
    SOC_soft = SOC_soft + SOC_penalty(7)*(NEAR_SOC_max(8)  > table_x1 & table_x1 > NEAR_SOC_max(7));
    SOC_soft = SOC_soft + SOC_penalty(8)*(NEAR_SOC_max(9)  > table_x1 & table_x1 > NEAR_SOC_max(8));
    SOC_soft = SOC_soft + SOC_penalty(9)*(NEAR_SOC_max(10) > table_x1 & table_x1 > NEAR_SOC_max(9));
    SOC_soft = SOC_soft + SOC_penalty(10)*(NEAR_SOC_max(11)> table_x1 & table_x1 > NEAR_SOC_max(10));
    SOC_soft = SOC_soft + SOC_penalty(11)*(NEAR_SOC_max(12)> table_x1 & table_x1 > NEAR_SOC_max(11));
    SOC_soft = SOC_soft + SOC_penalty(12)*(NEAR_SOC_max(13)> table_x1 & table_x1 > NEAR_SOC_max(12));
    SOC_soft = SOC_soft + SOC_penalty(13)*(NEAR_SOC_max(14)> table_x1 & table_x1 > NEAR_SOC_max(13));
    SOC_soft = SOC_soft + SOC_penalty(14)*(NEAR_SOC_max(15)> table_x1 & table_x1 > NEAR_SOC_max(14));
    SOC_soft = SOC_soft + SOC_penalty(15)*(NEAR_SOC_max(16)> table_x1 & table_x1 > NEAR_SOC_max(15));
    SOC_soft = SOC_soft + SOC_penalty(16)*(NEAR_SOC_max(17)> table_x1 & table_x1 > NEAR_SOC_max(16));
    SOC_soft = SOC_soft + SOC_penalty(17)*(NEAR_SOC_max(18)> table_x1 & table_x1 > NEAR_SOC_max(17));
    SOC_soft = SOC_soft + SOC_penalty(18)*(NEAR_SOC_max(19)> table_x1 & table_x1 > NEAR_SOC_max(18));
    SOC_soft = SOC_soft + SOC_penalty(19)*(NEAR_SOC_max(20)> table_x1 & table_x1 > NEAR_SOC_max(19));
    SOC_soft = SOC_soft + SOC_penalty(20)*((MAX_SOC > table_x1) & (table_x1 > NEAR_SOC_max(20)));
    
    inst_fuel = repmat(inst_fuel,[1,1,1,1,1,x1_length]); % Add an extra dimension for the fuel table
    inst_fuel = permute(inst_fuel,[6 1 2 3 4 5]);
    
    table_L = inst_fuel + SOC_soft + 10*infeasible_SOC + 200*infeasible_We + 200*infeasible_Tm + 200*infeasible_Wm + 200*infeasible_Te + 200*infeasible_Pbatt;   %[x2]x[u1]x[u2]x[u3]
    
    savename = ['Transitional Cost = ',num2str(t),' Table.mat'];
    save(savename,'table_x1','table_L');
    
    if RUN_TYPE == 0  % For DP only runs
        complete = 100  - (time_cyc - t)/(time_cyc)*100;
        clc
        fprintf('__________________________________________________\n\n')
        fprintf('Percent Complete of Dynamic Simulation = ')
        fprintf(num2str(complete))
        fprintf('\n')
        fprintf('__________________________________________________\n\n')
    end 
end
cd ..  % Come out of the folder
%%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%-------------------------Dynamic Programming-----------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Define Parameters
BETA = 91000;
Desired_SOC = 0.55; % When you do the optimization - Extract solutions from the middle

folder = [cyc_name, ' TABLES'];
cd(folder);      % Go into the tables that were previously made

J_STAR = repmat(BETA*(x1_grid - Desired_SOC).^2, [1, x2_length,x3_length]); % terminal state penalty, [x1]x[x2]

for t = time_cyc:-1:1
    loadfile_name = ['Transitional Cost = ',num2str(t),' TABLE.mat'];
    load(loadfile_name); 
    SOC_State_Penalty = single(zeros(x1_length,x2_length,x3_length,u1_length,u2_length,u3_length));
    
    for x2 = 1:x2_length
        for x3 = 1:x3_length
            for u1 = 1:u1_length
                for u2 = 1:u2_length
                    for u3 = 1:u3_length
                        
                        % Next Gear
                        if u2 == 1 && x3 == 1 || u2 == 3 && x3 == 6
                            u2_c = 0; % Cannot Shift
                        else
                            u2_c = u2_grid(u2);
                        end
                        x3_n = x3 + u2_c;
                        
                        % Next Engine State
                        if x2 == 1   % Engine off
                            u3_temp = [0 1];
                            u3_c = u3_temp(u3); % It would just turn it on?
                        else                              % How was it turning off the engine?
                            u3_temp = [-1 0];
                            u3_c = u3_temp(u3);
                        end
                        x2_n = x2 + u3_c;
                        
                        F  = griddedInterpolant(x1_grid,J_STAR(:,x2_n,x3_n),'linear');  % Penalizing where they land! - Just penalizing the SOC
                        SOC_State_Penalty(:,x2,x3,u1,u2,u3) = F(table_x1(:,x2,x3,u1,u2,u3));
                    end
                end
            end
        end
    end
    
    J_temp = table_L + SOC_State_Penalty;
    
    for x1 = 1:x1_length
        for x2 = 1:x2_length
            for x3 = 1:x3_length
                S = squeeze(J_temp(x1,x2,x3,:,:,:));
                [minS,idx] = min(S(:));
                [u1,u2,u3] = ind2sub(size(S),idx);
                opt_id_u1(x1,x2,x3) = u1;
                opt_id_u2(x1,x2,x3) = u2;
                opt_id_u3(x1,x2,x3) = u3;
                
                % Define the new optimum value
                opt_value(x1,x2,x3) = J_temp(x1,x2,x3,u1,u2,u3);
            end
        end
    end
    J_STAR = opt_value;   % Next terminal cost!
    
    savename=['Cost & Control = ',num2str(t),' TABLE.mat'];
    save(savename,'J_STAR','opt_id_u1','opt_id_u2','opt_id_u3');
    
    if RUN_TYPE == 0
        complete = (time_cyc - t)/(time_cyc)*100;
        clc
        fprintf('__________________________________________________\n\n')
        fprintf('Percent Complete of Dynamic Programming = ')
        fprintf(num2str(complete))
        fprintf('\n')
        fprintf('__________________________________________________\n\n')
    end
end
cd ..

%%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%--------------------------Simulate Final Runs----------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
if RUN_TYPE == 0
    clc
    min_t  = toc/60;
    fprintf('__________________________________________________________________\n\n')
    fprintf('Total Time (min) for Dynamics and Dynamic Programming = ')
    fprintf(num2str(min_t))
    fprintf('\n')
    fprintf('__________________________________________________________________\n\n')
end

% Pre-Allocate for Speed
fail_inner_SOC = zeros(1,time_cyc);
fail_inner_Te = zeros(1,time_cyc);
fail_inner_We = zeros(1,time_cyc);
fail_inner_Tm = zeros(1,time_cyc);
fail_inner_Wm = zeros(1,time_cyc);
fail_inner_Pbatt = zeros(1,time_cyc);

if RUN_TYPE == 0
    fprintf('-------------------------------------------------\n');
    fprintf('Simulating Final Run! ')
    fprintf('\n')
    fprintf('-------------------------------------------------\n');
end

folder = [cyc_name, ' TABLES'];
cd(folder);

% Pre-Allocate for Speed
SOC_final = zeros(1,time_cyc);
inst_fuel = zeros(1,time_cyc);
W_mot = zeros(1,time_cyc);
T_mot = zeros(1,time_cyc);
W_eng = zeros(1,time_cyc);
T_eng = zeros(1,time_cyc);
ENG_state = zeros(1,time_cyc);
GEAR = zeros(1,time_cyc);
U3_save = zeros(1,time_cyc);
U1_sim = zeros(1,time_cyc);
U2_sim = zeros(1,time_cyc);
U3_sim = zeros(1,time_cyc);
GEAR_save = zeros(1,time_cyc);
J = zeros(1,time_cyc);
Pbatt_sim = zeros(1,time_cyc);

% Define Initial Conditions
ENG_state_id = 1;               % Engine Off
ENG_state_c = x2_grid(ENG_state_id);
SOC_c = 0.55;
GEAR_id = 1;                    % Start in First Gear

for t = 1:1:time_cyc
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    %------------- Load & Determine all Control Signals --------------%
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    load(['Cost & Control = ',num2str(t),' TABLE.mat']);
    
    % Engine Torque Control
    id_lookup_u1 = interp1(x1_grid,opt_id_u1(:,ENG_state_id,GEAR_id),SOC_c,'nearest');
    
    % Shifting Control
    id_lookup_u2 = interp1(x1_grid,opt_id_u2(:,ENG_state_id,GEAR_id),SOC_c,'nearest');
    u2_c = u2_grid(id_lookup_u2);
    
    % Engine Control
    id_lookup_u3 =  interp1(x1_grid,opt_id_u3(:,ENG_state_id,GEAR_id),SOC_c,'nearest');
    
    J(t) =  interp1(x1_grid,J_STAR(:,ENG_state_id,GEAR_id),SOC_c,'nearest');
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    
    % Shifting Control Check
    if (GEAR_id == 1 && u2_c == u2_grid(1)) || (GEAR_id == 6 && u2_c == u2_grid(u2_length))
        u2_c = 0;             % Do nothing -Maybe remove this later and penalize or something
    end
    
    % Identify Engine Control Vector - Depends on engine state
    if ENG_state_c == 0;
        u3_temp = [0 1];
        u3_c = u3_temp(id_lookup_u3); % It would just turn it on?
    else                              % How was it turning off the engine?
        u3_temp = [-1 0];
        u3_c = u3_temp(id_lookup_u3);
    end
    
    % Update x2 and define new eng id 
    ENG_state_n = ENG_state_c + u3_c;
    
    % Update x3 and define new gear id
    New_Gear_Index = GEAR_id + u2_c;
    GEAR_n = x3_grid(New_Gear_Index);
    
    if ENG_state_n == 0;
        Te_c = 0;
        We_c = 0;                              % [rad/sec]
        Wm_c = Ww(t)*FD*G;              % [rad/sec]
    else
        We_c = Ww(t)*FD*GEAR_n;                % [rad/sec]
        Taux = Paux/We_c;
        Te_c = u1_grid(id_lookup_u1) - Taux;
        Wm_c = Ww(t)*FD*G;              % [rad/sec]
    end
    
    if Pd(t) < 0 && ENG_state_n == 1
        We_c = optimal_eng_spd;
        Taux = Paux/We_c;
        Te_c = Taux;
        Tm_c = Tw(t)/(FD*G);  % Engine is disconnected
    else
        Tm_c = Tw(t)/(FD*G) - Te_c*GEAR_n/G;  % [1]x[1]
    end
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    %----- Check prime movers & Saturate For Fuel & Eff. Tables ------%
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    %                         ~Engine Speed~
    if We_c > W_eng_max
        We_fuel = W_eng_max;
        Fail_We = 1;
    elseif We_c < W_eng_min
        We_fuel = W_eng_min;
        if We_c < 800*rpm2rads &&  ENG_state_n == 1 % Engine is on
            Fail_We = 1;
        else
            Fail_We = 0;  % Don't fail it if the engine is off
        end
    else
        We_fuel = We_c;
        Fail_We = 0;
    end
    
    Te_max = interp1(eng_map_spd,eng_max_trq,We_fuel);
    %                         ~Engine Torque~
    if Te_c > Te_max
        Te_fuel = Te_max;
        Fail_Te = 1;
    elseif Te_c < Te_min(1)
        Te_fuel = Te_min(1);
        Fail_Te = 0;
    else
        Te_fuel = Te_c;
        Fail_Te = 0;
    end
    
    if ENG_state_n == 1;
        fuel = interp2(eng_consum_trq,eng_consum_spd,eng_consum_fuel,Te_fuel,We_fuel,'linear')*dt;
    else
        fuel = 0;
    end
    
    %                          ~Motor Speed~
    if Wm_c > Wm_max
        Wm_eff = Wm_max;
        Fail_Wm = 1;
    elseif Wm_c < Wm_min
        Wm_eff = Wm_min;
        Fail_Wm = 1;
    else
        Wm_eff = Wm_c;
        Fail_Wm = 0;
    end
    
    Tm_max = interp1(m_map_spd,m_max_trq,Wm_eff);
    %                          ~Motor Torque~
    if Tm_c > Tm_max
        Tm_eff = Tm_max;
        Fail_Tm = 1;
    elseif Tm_c < -Tm_max
        Tm_eff = -Tm_max;
        Fail_Tm = 0;    % Can use the brake
    else
        Tm_eff = Tm_c;
        Fail_Tm = 0;
    end
    eff_m = interp2(m_map_trq, m_map_spd, m_eff_map, Tm_eff, abs(Wm_eff)); % Assume eff. table is symetric
    eff_m(isnan(eff_m)) = 0.2;
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    %------------- Calculate New SOC & Check New SOC------------------%
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    
    if Tm_c <= 0
        Pbat = (Wm_c*Tm_c)*(eff_m);
    else
        Pbat = (Wm_c*Tm_c)/(eff_m);
    end
    
    if ENG_state_n == 0;
        Pbat_terminal = Pbat + Paux;   % If Pbat is positive, we are adding to the system and taking away from the battery SOC
        % A positive Paux takes away from the SOC
    else
        Pbat_terminal = Pbat;
    end
    
    if Pbat_terminal >= 0
        Pbatt_max = interp1(ess_soc, ess_max_pwr_dis,SOC_c);
        if Pbat_terminal > Pbatt_max
            FAIL_Pbatt = 1;
            Pbat_terminal = Pbatt_max;  % Saturate it 
        else
            FAIL_Pbatt = 0;
        end
%         Pbat(Pbat > Pbatt_max) = Pbatt_max;
        rint_c = interp1(ess_soc,ess_r_dis,SOC_c);
    else
        Pbatt_min = -interp1(ess_soc, ess_max_pwr_chg,SOC_c);
        FAIL_Pbatt = 0;      % Brakes will handle the rest of the power
        Pbat_terminal(Pbat_terminal < Pbatt_min) = Pbatt_min;
        rint_c = interp1(ess_soc,ess_r_chg,SOC_c);
    end
    Voc_c = interp1(ess_soc,ess_voc,SOC_c);
    i(t) = (Voc_c-(Voc_c^2-4*Pbat_terminal*rint_c)^(1/2))/(2*rint_c); % Picked the smaller current (-)
    SOC_n = SOC_c - i(t)/(ess_cap_ah*3600)*dt;
%   SOC_n = SOC_c -(Voc_c-(Voc_c^2-4*Pbat_terminal*rint_c)^(1/2))/(2*rint_c*ess_cap_ah*3600)*dt;
    C_rate(t) = i(t)/ess_cap_ah;
    
    % Determine Braking Torque - With updated Pbat
    if Wm_c ~= 0
        if ENG_state_n == 0;  % Don't suppply the auxilary power though the motor
            if Tm_c <= 0 
                Tm_actual =  (Pbat_terminal - Paux)/(Wm_c*eff_m); 
            else
                Tm_actual = (Pbat_terminal - Paux)/(Wm_c)*eff_m;
            end
        else
            if Tm_c <= 0
                Tm_actual =  (Pbat_terminal)/(Wm_c*eff_m);
            else
                Tm_actual = (Pbat_terminal)/(Wm_c)*eff_m;
            end
        end
    else
        Tm_actual = 0;
    end
    T_brake = Tm_c -  Tm_actual;
    
    % Check new SOC
    if SOC_n > MAX_SOC
        SOC_n = MAX_SOC;
        Fail_SOC = 1;
    elseif SOC_n < MIN_SOC
        SOC_n = MIN_SOC;
        Fail_SOC = 1;
    else
        Fail_SOC = 0;
    end
    
    % Save the opperational feasiblilty results
    fail_inner_SOC(t) = Fail_SOC;
    fail_inner_Te(t) = Fail_Te;
    fail_inner_We(t) = Fail_We;
    fail_inner_Tm(t) = Fail_Tm;
    fail_inner_Wm(t) = Fail_Wm;
    fail_inner_Pbatt(t) = FAIL_Pbatt;
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    %---------------------- Update the States-------------------------%
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    % Update x1 for next time step
    SOC_c = SOC_n;
    
    % Update x2 for next time step
    ENG_state_c = ENG_state_n;
    ENG_state_id = find(ENG_state_n == x2_grid);
    
    % Update x3 for next time step
    GEAR_id = find(GEAR_n == x3_grid);
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    
    % Save Simulation Variables
    SOC_final(t) = SOC_c;
    ENG_state(t) = ENG_state_c;
    GEAR(t) = GEAR_id;
    GEAR_save(t) = GEAR_n; % Actual gear ratio
    W_eng(t) = We_c;
    T_eng(t) = Te_c;
    W_mot(t) = Wm_c;
    T_mot(t) = Tm_c;
    inst_fuel(t) = fuel;
    U3_save(t) = u3_c;
    Pbatt_sim(t) = Pbat_terminal;
    eff_m_sim(t) = eff_m;
    
    T_brake_sim(t) = T_brake;
    Tm_actual_sim(t) = Tm_actual;
    
    % Plot the control signal indicies
    U1_sim(t)  = id_lookup_u1;
    U2_sim(t)  = id_lookup_u2;
    U3_sim(t)  = id_lookup_u3;
end

T_eng(T_eng<0) = 0;
dSOC = SOC_final(end) - SOC_final(1);
total_fuel_gram = sum(inst_fuel);
total_distance_mile = sum(cyc_mph(:,2))/3600;
MPG = total_distance_mile/(total_fuel_gram/1000/gasoline_density*liter2gallon);
cd ..    % Come out of folder

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%-----------------------Check Outer Feasibility---------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%SOC
if any(any(fail_inner_SOC))
    fail_outer_SOC = 1;
else
    fail_outer_SOC = 0;
end
%We
if any(any(fail_inner_We))
    fail_outer_We = 1;
else
    fail_outer_We = 0;
end
%Te
if any(any(fail_inner_Te))
    fail_outer_Te = 1;
else
    fail_outer_Te = 0;
end
%Wm
if any(any(fail_inner_Wm))
    fail_outer_Wm = 1;
else
    fail_outer_Wm = 0;
end
%Tm
if any(any(fail_inner_Tm))
    fail_outer_Tm = 1;
else
    fail_outer_Tm = 0;
end
%Pbatt
if any(any(fail_inner_Pbatt))
    fail_outer_Pbatt = 1;
else
    fail_outer_Pbatt = 0;
end
FAIL = ((fail_outer_Tm + fail_outer_Wm + fail_outer_Te + fail_outer_We + fail_outer_SOC + fail_outer_Pbatt) ~= 0);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%-----------------------Final Plots ect.----------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%


if RUN_TYPE == 0
    cd('Plots')
    Main_Plot;
    Engine_Plot;
    Motor_Plot;
%     Cost_Plot;
%     Torque_Check;
%     Speed_Check;
    Battery_Plot;
%     Brake_Plot;
    cd ..
    MPG
    FAIL
end
%%

