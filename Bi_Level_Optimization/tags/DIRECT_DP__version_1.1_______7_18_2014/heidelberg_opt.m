close all;
clear all;
clc;
% Fixed engine state in final sinulation - 7/14 -Huck

% Load the Engine and Motor Data
cd('Components');
Engine_2rz_0410; 
Motor_int;
Battery_int;
cd ..  % Come up one folder

% Define Vehicle Parameters
Vehicle_Parameters;

%load CYC_HWFET; cyc_name = 'HWFET';
%load CYC_UDDS; cyc_name = 'UDDS';
%load CYC_US06; cyc_name = 'US06';
%load SHORT_CYC_HWFET; cyc_name = 'SHORT_CYC_HWFET';
load RAMP; cyc_name = 'RAMP';
%load RAMP_slow; cyc_name = 'RAMP_slow';

Manipulate_Drive_Cycle;

% Dynamic Programming

%Define State Grids
x1_grid = [0.4:0.05:0.8]';       % SOC 
x1_length = length(x1_grid);

x2_grid = [0 1];   % Engine off (0)   &   Engine on (1)
x2_length = length(x2_grid);

x3_grid = [4.484 2.872 1.842 1.414 1.000 0.742];  % [1st 2nd...]             % Gear Level
x3_length = length(x3_grid);

% Define Control Grids
u1_grid = [0:10:140]';          % Engine Torque [N-m]
u1_length = length(u1_grid);  

u2_grid = [-1 0 1];             % Shift down, nothing, up
u2_length = length(u2_grid);

u3_length = 2;

% Define Some soft cnstraints stuff
SOC_penalty = linspace(0.1,10,20);
NEAR_SOC_min = MIN_SOC + fliplr(linspace(0.001,0.02,20));
NEAR_SOC_max = MAX_SOC - fliplr(linspace(0.001,0.02,20));

% Design the Design Variables Range
FD_range = linspace(2,4,14);
G_range = linspace(2,4,14);
size_sim = length(FD_range) + length(G_range);
hour = 0;

tic
for ii = 1:length(FD_range)
    for rr = 1:length(G_range)
        FD = FD_range(ii);
        G = G_range(rr);
        
        %% Make some new folders - Should I Clear Them?
        tables = [cyc_name, ' TABLES'];
        mkdir(tables)
        cd(tables); % path into the folder
        
        for t = 1:time_cyc
            table_x1 = single(zeros(x1_length,x2_length,x3_length,u1_length,u2_length,u3_length));   % [x1]x[x2]x[x3]x[u1]x[u2]x[u3]
            table_L = single(zeros(x1_length,x2_length,x3_length,u1_length,u2_length,u3_length));    % [x1]x[x2]x[x3]x[u1]x[u2]x[u3]
            Tm_max = single(zeros(x2_length,x3_length,u1_length,u2_length,u3_length));               % [x2]x[u1]x[u2]
            Tm_save = single(zeros(x2_length,x3_length,u1_length,u2_length,u3_length));              % [x2]x[u1]x[u2]
            Wm_save = single(zeros(x2_length,x3_length,u1_length,u2_length,u3_length));              % [x2]x[u1]x[u2]
            We_save = single(zeros(x2_length,x3_length,u1_length,u2_length,u3_length));              % [x2]x[u1]x[u2]
            table_x1 = single(zeros(x1_length,x2_length,x3_length,u1_length,u2_length,u3_length));   % [x2]x[u1]x[u2]
            inst_fuel = single(zeros(x2_length,x3_length,u1_length,u2_length,u3_length));            % [x2]x[u1]x[u2]
            infeasible_Te = single(zeros(x2_length,x3_length,u1_length,u2_length,u3_length));        % [x2]x[u1]x[u2]
            for x3 = 1:x3_length           % Go through all of the gears
                x3_c = x3_grid(x3);
                
                for u2 = 1:u2_length       % Shift down, don't shift, and shift up
                    u2_c = u2_grid(u2);
                    
                    
                    if (x3_c == x3_grid(1) && u2_c == u2_grid(1)) || (x3_c == x3_grid(x3_length) && u2_c == u2_grid(u2_length))
                        u2_c = 0;          % Do nothing - Cannot shift!!
                    end
                    
                    if u2 == 1 || u2 == 3  % Shift penalty
                        Shift_Penalty = 0;
                    else
                        Shift_Penalty = 0;
                    end
                    
                    % Update x3 and define new gear id
                    New_Gear_Index = x3 + u2_c;
                    x3_n = x3_grid(New_Gear_Index);
                    
                    for x2 = 1:x2_length             % Current Engine State
                        ENG_state_c = x2_grid(x2);
                        
                        for u3 = 1:u3_length         % Engine Control
                            
                            if  ENG_state_c == 0 && u3 == 2   % Engine will be turned on - Done later in the code..this is the same thing though.
                                Eng_Penalty = 0;
                            else
                                Eng_Penalty = 0;
                            end
                            
                            if ENG_state_c == 0 && u3 == 1 || ENG_state_c == 1 && u3 == 1    % Next Engine State is off
                                ENG_state_n = 0;
                                We_c = 0;                            % [rad/sec]
                                Wm_c = Ww(t)*FD*x3_n*G;              % [rad/sec]
                                Te_c = zeros(u1_length,1);           % Turn the engine off
                            else                                                             % Engine is on
                                ENG_state_n = 1;
                                Te_c = u1_grid;                      % Engine Control, [u1]x[1]
                                We_c = Ww(t)*FD*x3_n;                % [rad/sec]
                                Wm_c = Ww(t)*FD*x3_n*G;              % [rad/sec]
                            end
                            
                            if Pd(t) < 0 && ENG_state_n == 1 % Braking and Running Engine - Idle Engine
                                Te_c = zeros(u1_length,1);                         % Engine Does not add any torque!
                                We_c = 1000*rpm2rads;                              % Assumed idle speed ( rads/sec )
                            end
                            
                            Tm_c = Tw(t)/(x3_n*FD*G)*ones(size(Te_c)) - Te_c/G;  % [u1]x[1]
                            
                            % Check Motor
                            Tm_max(x2,x3,:,u2,u3) = interp1(m_map_spd,m_max_trq,abs(Wm_c))*ones(size(Tm_c));    % [x2]x[x3]x[u1]x[u2]x[u3]   - #Check to see if the motor map is symmetric
                            Tm_save(x2,x3,:,u2,u3) = Tm_c;                                                      % [x2]x[x3]x[u1]x[u2]x[u3]
                            Wm_save(x2,x3,:,u2,u3) = Wm_c*ones(size(u1_grid));                                  % [u1]x[1]  -  Check For Each Gear
                            
                            % Check Engine
                            We_save(x2,x3,:,u2,u3) = We_c*ones(size(u1_grid));                         % [u1]x[1]  - Check For Each Gear - Speed of the engine depends on the gear and engine state
                            Te_max = interp1(eng_map_spd,eng_max_trq,abs(We_c))*ones(size(Te_c));
                            infeasible_Te(x2,x3,:,u2,u3) = (Te_max < Te_c);                            % [u1]x[1]
                            
                            % Update x1
                            eff_m = interp2(m_map_trq, m_map_spd, m_eff_map, Tm_c, abs(Wm_c))';        % [u1]x[1]
                            eff_m(isnan(eff_m)) = 0.2;
                            Pm_c = (Wm_c*Tm_c).*(eff_m.^-sign(Wm_c*Tm_c));                             % If the torque is negative, we are charging the battery
                            Pbat = repmat(Pm_c,[1, x1_length]);                                        % [u1]x[1]
                            Pbat = permute(Pbat, [2 1]);
                            
                            if Pbat >= 0
                                Pbatt_max = repmat(interp1(ess_soc, ess_max_pwr_dis, x1_grid),[1,u1_length]);
                                Pbat(Pbatt_max < Pbat) = Pbatt_max(Pbatt_max < Pbat);
                                rint_c = repmat(interp1(ess_soc,ess_r_dis,x1_grid),[1,u1_length]);
                            else
                                Pbatt_min = repmat(interp1(ess_soc, ess_max_pwr_chg, x1_grid),[1,u1_length]);
                                Pbat(-Pbatt_min > Pbat) = -Pbatt_min(-Pbatt_min > Pbat);  % Pbatt_min hould be more negative
                                rint_c = repmat(interp1(ess_soc,ess_r_chg,x1_grid),[1,u1_length]);
                            end
                            voc_c = repmat(interp1(ess_soc,ess_voc,x1_grid), [1, u1_length]);
                            SOC_c_matrix = repmat(x1_grid,[1, u1_length]);
                            SOC_n = real((-(voc_c-(voc_c.^2-4*Pbat.*rint_c).^(1/2))./(2*rint_c)/(ess_cap_ah*3600))*dt + SOC_c_matrix);  %[u1]x[1]
                            Out_Of_Range = 10;  % Can move outside for loop later
                            fuel = (interp2(eng_consum_trq,eng_consum_spd,eng_consum_fuel,Te_c,We_c,'linear',Out_Of_Range)*dt)';                 %[1]x[u1]
                            table_x1(:,x2,x3,:,u2,u3) = SOC_n;                                      %[x2]x[x3]x[u1]x[u2]x[u3]- Different SOC's depending on the control
                            inst_fuel(x2,x3,:,u2,u3) = fuel + Shift_Penalty*ones(size(fuel)) + Eng_Penalty*ones(size(fuel));                      %[x2]x[x3]x[u1]x[u2]x[u3]- Different fuel's depending on the control
                            
                        end   % End of u3 ( Engine Control ) Loops
                    end       % End of x2 ( Engine State ) Loops
                end           % End of u2 ( Gear Control )
            end               % End of x3 (Gear State) loops
            
            % Check Motor
            infeasible_Tm = (Tm_save > Tm_max) | (Tm_save < -Tm_max);             % [x2]x[x3]x[u1]x[u2]x[u3]
            infeasible_Tm = repmat(infeasible_Tm,[1,1,1,1,1,x1_length]);
            infeasible_Tm = permute(infeasible_Tm,[6 1 2 3 4 5]);
            
            infeasible_Wm = (Wm_save > Wm_max) | (Wm_save < Wm_min);              % [x2]x[x3]x[u1]x[u2]x[u3]
            infeasible_Wm = repmat(infeasible_Wm,[1,1,1,1,1,x1_length]);
            infeasible_Wm = permute(infeasible_Wm,[6 1 2 3 4 5]);
            
            % Check Engine - Torque is defined, so it will be ok. If we had engine speed as a state, we could check the acceleration
            infeasible_We = (We_save < 800*rpm2rads) | (We_save > W_eng_max);    % [x2]x[x3]x[u1]x[u2]
            infeasible_We = repmat(infeasible_We,[1,1,1,1,1,x1_length]);
            infeasible_We = permute(infeasible_We,[6 1 2 3 4 5]);
            
            infeasible_Te = repmat(infeasible_Te,[1,1,1,1,1,x1_length]);
            infeasible_Te = permute(infeasible_Te,[6 1 2 3 4 5]);
            
            % Check SOC
            infeasible_SOC = (table_x1 < MIN_SOC) | (table_x1 > MAX_SOC);        % [x2]x[x3]x[u1]x[u2]x[u3]
            table_x1(table_x1 > MAX_SOC) = MAX_SOC;  % Check this as well!!
            table_x1(table_x1 < MIN_SOC) = MIN_SOC;
            
            SOC_soft = zeros(size(table_x1));                                    % [x2]x[x3]x[u1]x[u2]x[u3]
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
            
            % Add an extra dimension for the fuel table
            inst_fuel = repmat(inst_fuel,[1,1,1,1,1,x1_length]);
            inst_fuel = permute(inst_fuel,[6 1 2 3 4 5]);
            
            table_L = inst_fuel + SOC_soft + 10*infeasible_SOC + 10*infeasible_We + 10*infeasible_Tm + 10*infeasible_Wm + 10*infeasible_Te;   %[x2]x[u1]x[u2]x[u3]
            
            savename = ['Transitional Cost = ',num2str(t),' Table.mat'];
            save(savename,'table_x1','table_L');
            
            % Display Status
            complete = 100  - (time_cyc - t)/(time_cyc)*100;
            total_sim = (rr + ii - 2)/size_sim*100;
            clc
            fprintf('__________________________________________________\n\n')
            fprintf('Percent Complete of Dynamic Simulation = ')
            fprintf(num2str(complete))
            fprintf('\n')
            fprintf('__________________________________________________\n\n\n\n')
            fprintf('__________________________________________________\n')
            fprintf('Percent Complete of Total Simulation = ')
            fprintf(num2str(total_sim))
            fprintf('\n')
            fprintf('__________________________________________________\n\n')
        end
        cd ..  % Come out of the folder
        %%
        % Define Parameters
        BETA = 9100;
        Desired_SOC = 0.55; % When you do the optimization - Extract solutions from the middle
        
        folder = [cyc_name, ' TABLES'];
        cd(folder);      % Go into the tables that were previously made
        
        J_STAR = repmat(BETA*(x1_grid - Desired_SOC).^2, [1, x2_length,x3_length]); % terminal state penalty, [x1]x[x2]
        
        for t = time_cyc:-1:1
            loadfile_name = ['Transitional Cost = ',num2str(t),' TABLE.mat'];
            load(loadfile_name);
            
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
                        
                        clear S;
                    end
                end
            end
            
            J_STAR = opt_value;   % Next terminal cost!
            
            savename=['Cost & Control = ',num2str(t),' TABLE.mat'];
            save(savename,'J_STAR','opt_id_u1','opt_id_u2','opt_id_u3');
            
            % Display Status
            complete = (time_cyc - t)/(time_cyc)*100;
            clc
            fprintf('__________________________________________________\n\n')
            fprintf('Percent Complete of Dynamic Programming = ')
            fprintf(num2str(complete))
            fprintf('\n')
            fprintf('__________________________________________________\n\n')
        end
        cd ..
        %%
        hour  = hour + toc/3600;
        clc
        fprintf('__________________________________________________________________\n\n')
        fprintf('Total Time Elapsed (hours) = ')
        fprintf(num2str(hour))
        fprintf('\n')
        fprintf('__________________________________________________________________\n\n')
        
        u = 1;
        for SOC_c = 0.45:0.05:0.75
            
            folder = [cyc_name, ' TABLES'];
            cd(folder);
            
            % Pre allocate for speed
            SOC_final = zeros(1, time_cyc);
            inst_fuel = zeros(1,time_cyc);
            W_mot = zeros(1, time_cyc);
            T_mot = zeros(1, time_cyc);
            W_eng = zeros(1, time_cyc);
            T_eng = zeros(1,time_cyc);
            ENG_state = zeros(1,time_cyc);
            GEAR = zeros(1,time_cyc);
            U3_save = zeros(1,time_cyc);
            U1_sim = zeros(1,time_cyc);
            U2_sim = zeros(1,time_cyc);
            U3_sim = zeros(1,time_cyc);
            GEAR_save = zeros(1,time_cyc);
            J = zeros(1,time_cyc);
            
            % Define Initial Conditions
            ENG_state_id = 1;               % Engine Off
            ENG_state_c = x2_grid(ENG_state_id);
            
            GEAR_id = 1;                    % Start in First Gear
            
            for t = 1:1:time_cyc
                
                load(['Cost & Control = ',num2str(t),' TABLE.mat']);
                
                %_____________________Load all control signals________________________%
                
                % Engine Torque Control
                id_lookup_u1 = interp1(x1_grid,opt_id_u1(:,ENG_state_id,GEAR_id),SOC_c,'nearest');
                
                % Shifting Control
                id_lookup_u2 = interp1(x1_grid,opt_id_u2(:,ENG_state_id,GEAR_id),SOC_c,'nearest');
                u2_c = u2_grid(id_lookup_u2);
                
                % Engine Control
                id_lookup_u3 =  interp1(x1_grid,opt_id_u3(:,ENG_state_id,GEAR_id),SOC_c,'nearest');
                
                %_____________________________________________________________________%
                
                J(t) =  interp1(x1_grid,J_STAR(:,ENG_state_id,GEAR_id),SOC_c,'nearest');
                
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
                
                % Update x2 and define new eng id  - change to new for consitiancy
                ENG_state_n = ENG_state_c + u3_c;
                
                % Update x3 and define new gear id
                New_Gear_Index = GEAR_id + u2_c;
                GEAR_n = x3_grid(New_Gear_Index);
                
                if ENG_state_n == 0;
                    Te_c = 0;
                    We_c = 0;                              % [rad/sec]
                    Wm_c = Ww(t)*FD*GEAR_n*G;              % [rad/sec]
                else
                    Te_c = u1_grid(id_lookup_u1);
                    We_c = Ww(t)*FD*GEAR_n;                % [rad/sec]
                    Wm_c = Ww(t)*FD*GEAR_n*G;              % [rad/sec]
                end
                
                if Pd(t) < 0 && ENG_state_n == 1   % Stopped and the Engine is running
                    Te_c = 0;
                    We_c = 1000*rpm2rads;           % Assumed idle speed ( rads/sec )
                end
                
                Tm_c = Tw(t)/(GEAR_n*FD*G) - Te_c/G;  % [1]x[1]
                Tm_max = interp1(m_map_spd,m_max_trq,abs(Wm_c));
                Te_max = interp1(eng_map_spd,eng_max_trq,abs(We_c));
                %__________________________________________________________________%
                %-------------------- Check all prime movers ----------------------%
                
                %------Motor--------------Torque
                if Tm_c > Tm_max
                    Tm_c = Tm_max;
                    Fail_Tm = 1;
                elseif Tm_c < -Tm_max
                    Tm_c = -Tm_max;
                    Fail_Tm = 1;
                else
                    Fail_Tm = 0;
                end
                %------Motor--------------Speed
                if Wm_c > Wm_max
                    Wm_c = Wm_max;
                    Fail_Wm = 1;
                elseif Wm_c < -Wm_max
                    Wm_c = -Wm_max;
                    Fail_Wm = 1;
                else
                    Fail_Wm = 0;
                end
                
                %------Engine--------------Torque
                if Te_c > Te_max
                    Te_c = Te_max;
                    Fail_Te = 1;
                else
                    Fail_Te = 0;
                end
                %------Engine--------------Speed
                if We_c > W_eng_max
                    We_c = W_eng_max;
                    Fail_We = 1;
                elseif 0 > We_c  % Cannot do 800 RPM?
                    We_c = 0;
                    Fail_We = 1;
                else
                    Fail_We = 0;
                end
                %-------------------- Check all prime movers ----------------------%
                %__________________________________________________________________%
                
                fuel = interp2(eng_consum_trq,eng_consum_spd,eng_consum_fuel,Te_c,We_c,'linear')*dt;   %[1]x[1]
                
                eff_m = interp2(m_map_trq, m_map_spd, m_eff_map, Tm_c, abs(Wm_c)); % Check this, the motor efficiency is probably not symetric...could mess up Pbat
                eff_m(isnan(eff_m)) = 0.2; % to avoid nan
                Pm_c = (Wm_c*Tm_c)*(eff_m^-sign(Wm_c*Tm_c));
                Pbat = Pm_c;     % [1]x[1] - fixed 7/14
                
                if Pbat >= 0
                    Pbatt_max = interp1(ess_soc, ess_max_pwr_dis,SOC_c);
                    Pbat(Pbat > Pbatt_max) = Pbatt_max;
                    rint_c = interp1(ess_soc,ess_r_dis,SOC_c);
                else
                    Pbatt_min = interp1(ess_soc, ess_max_pwr_chg,SOC_c);
                    Pbat(Pbat < -Pbatt_min) = -Pbatt_min;
                    rint_c = interp1(ess_soc,ess_r_chg,SOC_c);
                end
                voc_c = interp1(ess_soc,ess_voc,SOC_c);
                SOC_n = ((-(voc_c-(voc_c^2-4*Pbat*rint_c)^(1/2))/(2*rint_c)/(ess_cap_ah*3600))*dt + SOC_c);
                
                %----------------------Check SOC----------------------------------%
                if SOC_n > MAX_SOC
                    SOC_n = MAX_SOC;
                    Fail_SOC = 1;
                elseif SOC_n < MIN_SOC
                    SOC_n = MIN_SOC;
                    Fail_SOC = 1;
                else
                    Fail_SOC = 0;
                end
                %----------------------Check SOC-----------------------------------%
                
                % Save the opperational feasiblilty results 
                fail_inner_SOC(u,t) = Fail_SOC;
                fail_inner_Te(u,t) = Fail_Te;
                fail_inner_We(u,t) = Fail_We;
                fail_inner_Tm(u,t) = Fail_Tm;
                fail_inner_Wm(u,t) = Fail_Wm;
                
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
                
                % Plot the control signal indicies
                U1_sim(t)  = id_lookup_u1;
                U2_sim(t)  = id_lookup_u2;
                U3_sim(t)  = id_lookup_u3;
                
                %__________________________Update States__________________________%
                % Update x1 for next time step
                SOC_c = SOC_n;
                
                % Update x2 for next time step
                ENG_state_c = ENG_state_n;
                ENG_state_id = find(ENG_state_n == x2_grid);
                
                % Update x3 for next time step
                GEAR_id = find(GEAR_n == x3_grid);
                
                %_________________________________________________________________%
            end
                       
            T_eng(T_eng<0) = 0;
            dSOC = SOC_final(end) - SOC_final(1);
            total_fuel_gram = sum(inst_fuel);
            total_distance_mile = sum(cyc_mph(:,2))/3600;
            mpg = total_distance_mile/(total_fuel_gram/1000/gasoline_density*liter2gallon);
            cd ..    % Come out of folder
            
            % SOC Correction
            MPG_save(u) = mpg;
            Delta_SOC(u) = dSOC;
            u = u+1;
            
        end
        %________________ Check all Outer Feasibilites____________________%
        %SOC
        if any(any(fail_inner_SOC))
            fail_outer_SOC(rr,ii) = 1;
        else
            fail_outer_SOC(rr,ii) = 0;
        end
        %We
        if any(any(fail_inner_We))
            fail_outer_We(rr,ii) = 1;
        else
            fail_outer_We(rr,ii) = 0;
        end
        %Te
        if any(any(fail_inner_Te))
            fail_outer_Te(rr,ii) = 1;
        else
            fail_outer_Te(rr,ii) = 0;
        end
        %Wm
        if any(any(fail_inner_Wm))
            fail_outer_Wm(rr,ii) = 1;
        else
            fail_outer_Wm(rr,ii) = 0;
        end
        %Tm
        if any(any(fail_inner_Tm))
            fail_outer_Tm(rr,ii) = 1;
        else
            fail_outer_Tm(rr,ii) = 0;
        end
        
        %________________ Check all Outer Feasibilites____________________%
        
        MPG_save = MPG_save';
        Delta_SOC = Delta_SOC';
        
        num = u -1;
        fit1 = polyfit(Delta_SOC(1:num), MPG_save(1:num), 1);
        SOC_Corrected_MPG = fit1(2);
        
        MPG(rr,ii) = SOC_Corrected_MPG;
        G_save(rr,ii) = G;
        FD_save(rr,ii) = FD;              
        
    end
end

Gear_Plot;

Main_Plot;

% %%
% figure(2);clf;
% [AX,H1,H2] = plotyy(cyc_time,cyc_mph(:,2),cyc_time,SOC_final);
% set(get(AX(1),'Ylabel'),'String','Vehicle Speed (mph)','fontWeight','bold','fontSize',20)
% set(get(AX(2),'Ylabel'),'String','SOC','fontWeight','bold','fontSize',20)
% set(H1,'LineWidth',2);
% set(H2,'LineWidth',2);
% set(AX(2),'fontWeight','bold','fontSize',20)
% xlabel('time (sec)','fontWeight','bold','fontSize',20);
% set(gca,'fontSize',20,'fontWeight','bold'),grid
%        
% %%
% figure(3);clf;
% [AX,H1,H2] = plotyy(cyc_time,cyc_mph(:,2),cyc_time,GEAR);
% set(get(AX(1),'Ylabel'),'String','Vehicle Speed (mph)','fontWeight','bold','fontSize',20)
% set(get(AX(2),'Ylabel'),'String','Optimal Gear Level','fontWeight','bold','fontSize',20)
% set(H1,'LineWidth',2);
% set(H2,'LineWidth',2);
% set(AX(2),'fontWeight','bold','fontSize',20)
% xlabel('time (sec)','fontWeight','bold','fontSize',20);
% set(gca,'fontSize',20,'fontWeight','bold'),grid
% 
% %%
% figure(4);clf;
% [AX,H1,H2] = plotyy(cyc_time,cyc_mph(:,2),cyc_time,inst_fuel);
% set(get(AX(1),'Ylabel'),'String','Vehicle Speed (mph)','fontWeight','bold','fontSize',20)
% set(get(AX(2),'Ylabel'),'String','Instantaneous Fuel (g/s)','fontWeight','bold','fontSize',20)
% set(H1,'LineWidth',2);
% set(H2,'LineWidth',2);
% set(AX(2),'fontWeight','bold','fontSize',20)
% xlabel('time (sec)','fontWeight','bold','fontSize',20);
% set(gca,'fontSize',20,'fontWeight','bold'),grid
% 
% %%
% figure(5);clf;
% plot(cyc_time,W_eng*rads2rpm,'color','g','LineWidth',2);
% hold on 
% plot(cyc_time,W_mot*rads2rpm,'LineWidth',2);grid on;
% legend('Engine Speed','Motor Speed')
% ylabel('Speed (rpm)','fontWeight','bold','fontSize',20)
% xlabel('time (sec)','fontWeight','bold','fontSize',20);
% set(gca,'fontSize',20,'fontWeight','bold')
% % title('HWFET cycle','fontWeight','bold','fontSize',20)
% 
% %%
% figure(6);clf;
% plot(cyc_time,T_eng,'color','g','LineWidth',2);
% hold on 
% plot(cyc_time,T_mot,'LineWidth',2);
% 
% legend('Engine Torque (Nm)','Motor Torque (Nm)')
% ylabel('Speed (rpm)','fontWeight','bold','fontSize',20)
% xlabel('time (sec)','fontWeight','bold','fontSize',20);
% set(gca,'fontSize',20,'fontWeight','bold')
% % title('HWFET cycle','fontWeight','bold','fontSize',20)
% 
% %% Check
% figure(7);clf;
% [AX,H1,H2] = plotyy(cyc_time,Pd,cyc_time,W_eng*rads2rpm);
% set(get(AX(1),'Ylabel'),'String','Pd (W)','fontWeight','bold','fontSize',20)
% set(get(AX(2),'Ylabel'),'String','Engine Speed (RPM)','fontWeight','bold','fontSize',20)
% set(H1,'LineWidth',2);
% set(H2,'LineWidth',2);
% set(AX(2),'fontWeight','bold','fontSize',20)
% xlabel('time (sec)','fontWeight','bold','fontSize',20);
% set(gca,'fontSize',20,'fontWeight','bold'),grid
%     
% 
% figure(8);clf;
% [AX,H1,H2] = plotyy(cyc_time,GEAR,cyc_time,W_eng*rads2rpm);
% set(get(AX(1),'Ylabel'),'String','Gear Level','fontWeight','bold','fontSize',20)
% set(get(AX(2),'Ylabel'),'String','Engine Speed (RPM)','fontWeight','bold','fontSize',20)
% set(H1,'LineWidth',2);
% set(H2,'LineWidth',2);
% set(AX(2),'fontWeight','bold','fontSize',20)
% xlabel('time (sec)','fontWeight','bold','fontSize',20);
% set(gca,'fontSize',20,'fontWeight','bold'),grid   
%  
% %%  
% figure(9);clf;
% [AX,H1,H2] = plotyy(cyc_time,ENG_state,cyc_time,W_eng*rads2rpm);
% set(get(AX(1),'Ylabel'),'String','Engine State','fontWeight','bold','fontSize',20)
% set(get(AX(2),'Ylabel'),'String','Engine Speed (RPM)','fontWeight','bold','fontSize',20)
% set(H1,'LineWidth',2);
% set(H2,'LineWidth',2);
% set(AX(2),'fontWeight','bold','fontSize',20)
% xlabel('time (sec)','fontWeight','bold','fontSize',20);
% set(gca,'fontSize',20,'fontWeight','bold'),grid   
%     
% %%   
% figure(10);clf;
% [AX,H1,H2] = plotyy(cyc_time,ENG_state,cyc_time,U3_save);
% set(get(AX(1),'Ylabel'),'String','Engine State','fontWeight','bold','fontSize',20)
% set(get(AX(2),'Ylabel'),'String','Engine Control','fontWeight','bold','fontSize',20)
% set(H1,'LineWidth',2);
% set(H2,'LineWidth',2);
% set(AX(2),'fontWeight','bold','fontSize',20)
% xlabel('time (sec)','fontWeight','bold','fontSize',20);
% set(gca,'fontSize',20,'fontWeight','bold'),grid   
%        
%     
% %%
% figure(11);clf;
% [AX,H1,H2] = plotyy(cyc_time,ENG_state,cyc_time,Pd);
% set(get(AX(1),'Ylabel'),'String','Engine State','fontWeight','bold','fontSize',20)
% set(get(AX(2),'Ylabel'),'String','Power Demand (W)','fontWeight','bold','fontSize',20)
% set(H1,'LineWidth',2);
% set(H2,'LineWidth',2);
% set(AX(2),'fontWeight','bold','fontSize',20)
% xlabel('time (sec)','fontWeight','bold','fontSize',20);
% set(gca,'fontSize',20,'fontWeight','bold'),grid      
%   
% %%
% figure(12); clf;
% subplot(3,1,1);
% GRAPH =  plot(cyc_time,U1_sim,'LineWidth',2); grid on;
% set(GRAPH,'marker','x', 'markersize', 6, 'markerf', 'b','linewidth',4);
% legend('From 1 to 15')
% ylabel('Engine Control - U1','fontWeight','bold','fontSize',12)
% xlabel('time (sec)','fontWeight','bold','fontSize',12);
% set(gca,'fontSize',12,'fontWeight','bold');   
% 
% subplot(3,1,2);
% [AX,L1,L2] = plotyy(cyc_time, U2_sim, cyc_time,GEAR);
% set(L1,'marker','x', 'markersize', 6, 'markerf','b','linewidth',4);
% set(L2,'marker','x', 'markersize', 6, 'markerf','g','linewidth',4);
% set(AX, 'fontsize', 12,'fontweight','bold');
% ylabel('Gear Control - U2', 'parent', AX(1));
% ylabel('Gear', 'parent', AX(2));grid
% xlabel('time (sec)');
% legend('[1 = DOWN] [2 = NOTHING] [3 = UP]','Gear Level');
% 
% subplot(3,1,3);
% [AX,L1,L2] = plotyy(cyc_time, U3_save, cyc_time,ENG_state);
% set(L1,'marker','x', 'markersize', 6, 'markerf', 'b','linewidth',4);
% set(L2,'marker','x', 'markersize', 2, 'markerf', 'g','linewidth',2);
% set(AX, 'fontsize', 12,'fontweight','bold');
% ylabel('Engine Control - U3', 'parent', AX(1));
% ylabel('Engine State (on/off)', 'parent', AX(2));grid
% xlabel('time (sec)');
% legend('[0 1] & [-1 0]','OFF  & ON')
% 
% %%
% Te = T_eng.*GEAR_save*FD;
% Tm = T_mot.*GEAR_save*FD*G;
% T_total = Te + Tm;
% 
% figure(13);clf;
% plot(cyc_time,Te,'color','g','LineWidth',2);
% hold on 
% plot(cyc_time,Tm,'LineWidth',3);
% hold on
% plot(cyc_time,Tw,'color','y','LineWidth',6);
% hold on
% plot(cyc_time,T_total,'r--','LineWidth',2)
% legend('Engine (at road)','Motor (at road)','T required','Motor + Engine')
% ylabel('Torque (Nm)','fontWeight','bold','fontSize',20);grid
% xlabel('time (sec)','fontWeight','bold','fontSize',20);
% set(gca,'fontSize',20,'fontWeight','bold')
% hold off    
%     
%     
% %%  Motor is always connected to vehicle...
% figure(14);clf;
% 
% We = W_eng./(GEAR_save*FD);
% Wm = W_mot./(GEAR_save*FD*G);
% 
% plot(cyc_time,Ww,'color','r','LineWidth',8);
% hold on 
% plot(cyc_time,Wm,'k--','LineWidth',4);
% hold on
% plot(cyc_time,We,'color','g','LineWidth',2);
% hold on
% legend('Wheel','Motor (at road)','Engine (at road)');
% ylabel('Speed (rad/s)','fontWeight','bold','fontSize',20)
% xlabel('time (sec)','fontWeight','bold','fontSize',20);
% set(gca,'fontSize',20,'fontWeight','bold'),grid
% hold off
% 
% %%
% figure(15);clf;
% [AX,H1,H2] = plotyy(cyc_time,Wm,cyc_time,GEAR);
% set(get(AX(1),'Ylabel'),'String','Motor Speed at road (rad/sec)','fontWeight','bold','fontSize',20)
% set(get(AX(2),'Ylabel'),'String','Gear Level','fontWeight','bold','fontSize',20)
% set(H1,'LineWidth',2);
% set(H2,'LineWidth',2);
% set(AX(2),'fontWeight','bold','fontSize',20)
% xlabel('time (sec)','fontWeight','bold','fontSize',20);
% set(gca,'fontSize',20,'fontWeight','bold'),grid
% 
% %%
% figure(16); clf;
% [C,h] = contour(eng_consum_spd*rads2rpm, eng_consum_trq, eng_bsfc', [216 200:5:250, 275, 300:100:800]);
% clabel(C,h, 'fontsize', 8);
% axis([0 5500 0 150]);
% hold on;
% plot(eng_map_spd*rads2rpm, eng_max_trq, 'r', 'linewidth', 4);
% % plot(We_optbsfc*rads2rpm, Te_optbsfc, 'ro-', 'markersize', 3, 'markerf', 'r');
% hold on 
% plot(W_eng*rads2rpm,T_eng, 'ko', 'markersize', 12, 'markerf', 'g')
% xlabel('Engine Speed (RPM)');
% ylabel('Torque (Nm)');
% set(gca,'FontSize',20,'fontWeight','bold')
% set(findall(gcf,'type','text'),'FontSize',25,'fontWeight','bold')
% legend('BSFC (g/kWh)','Maximum Engine Torque','Engine Opperating Points')    
%     
% %%    
% figure(17); clf;
% [C,h] = contour(m_map_spd*rads2rpm, m_map_trq, m_eff_map');
% 
% % axis([0 5500 0 150]);
% hold on; 
% plot(W_mot*rads2rpm,T_mot, 'ko', 'markersize', 8, 'markerf', 'g')
% hold on;
% plot(m_max_spd*rads2rpm, m_max_trq,  'ro', 'markersize', 12, 'markerf', 'r')
% hold on
% plot(m_max_spd*rads2rpm, m_max_gen_trq,  'ro', 'markersize', 12, 'markerf', 'r')
% xlabel('Motor Speed (RPM)');
% ylabel('Torque (Nm)');
% set(gca,'FontSize',20,'fontWeight','bold')
% set(findall(gcf,'type','text'),'FontSize',25,'fontWeight','bold')
% legend('Motor Efficiency','Motor Opperating Points','Maximum Torque'),grid 
% hold on
% clabel(C,h, 'fontsize', 25,'fontWeight','bold');    
% hold off    

%% 
% figure(19);clf;
% plot(cyc_time,J,'color','g','LineWidth',2);
% legend('Optimal Cost - With Penalty');
%     
% %% 
% figure(20); clf;
% subplot(3,1,1);
% GRAPH =  plot(cyc_time,U1_sim, 'kx-','LineWidth',2,'markersize', 3); grid on;
% set(GRAPH,'marker','x', 'markersize', n, 'markerf', 'b','linewidth',q);
% H = legend('Torque Control');set(H,'fontsize', h)
% ylabel({'Torque Control ID'},'fontWeight','bold','fontSize',y)
% set(gca,'fontSize',y,'fontWeight','bold');   
% set(gca,'XTickLabel',[]);
% xlim([0 t]);
% 
% 
% subplot(3,1,2);
% GRAPH =  plot(cyc_time,U2_sim, 'kx-','LineWidth',2,'markersize', 3); grid on;
% set(GRAPH,'marker','x', 'markersize', n, 'markerf', 'b','linewidth',q);
% H = legend('Gear Control');set(H,'fontsize', h)
% ylabel({'Gear Control ID'},'fontWeight','bold','fontSize',y)
% set(gca,'fontSize',y,'fontWeight','bold');   
% set(gca,'XTickLabel',[]);
% xlim([0 t]);
% 
% 
% subplot(3,1,3);
% GRAPH =  plot(cyc_time,U3_sim, 'kx-','LineWidth',2,'markersize', 3); grid on;
% set(GRAPH,'marker','x', 'markersize', n, 'markerf', 'b','linewidth',q);
% H = legend('Engine Control');set(H,'fontsize', h)
% ylabel({'Engine Control ID'},'fontWeight','bold','fontSize',y)
% set(gca,'fontSize',y,'fontWeight','bold');   
% set(gca,'XTickLabel',[]);
% xlim([0 t]);
% 
% %%
% figure(21); clf;
% subplot(3,1,1);
% GRAPH =  plot(cyc_time,SOC_final, 'kx-','LineWidth',2,'markersize', 3); grid on;
% set(GRAPH,'marker','x', 'markersize', n, 'markerf', 'b','linewidth',q);
% H = legend('SOC State');set(H,'fontsize', h)
% ylabel({'SOC'},'fontWeight','bold','fontSize',y)
% set(gca,'fontSize',y,'fontWeight','bold');   
% set(gca,'XTickLabel',[]);
% xlim([0 t]);
% 
% 
% subplot(3,1,2);
% GRAPH =  plot(cyc_time,GEAR, 'kx-','LineWidth',2,'markersize', 3); grid on;
% set(GRAPH,'marker','x', 'markersize', n, 'markerf', 'b','linewidth',q);
% H = legend('Gear State');set(H,'fontsize', h)
% ylabel({'Gear '},'fontWeight','bold','fontSize',y)
% set(gca,'fontSize',y,'fontWeight','bold');   
% set(gca,'XTickLabel',[]);
% xlim([0 t]);
% 
% 
% subplot(3,1,3);
% GRAPH =  plot(cyc_time,ENG_state, 'kx-','LineWidth',2,'markersize', 3); grid on;
% set(GRAPH,'marker','x', 'markersize', n, 'markerf', 'b','linewidth',q);
% H = legend('Engine State');set(H,'fontsize', h)
% ylabel({'Engine '},'fontWeight','bold','fontSize',y)
% set(gca,'fontSize',y,'fontWeight','bold');   
% set(gca,'XTickLabel',[]);
% xlim([0 t]);
%%
% figure(22);clf
% plot(Delta_SOC,MPG_save,'g.','markersize',35)
% hold on
% plot(0,fit1(2),'r.','markersize',35)
% title('SOC correction','fontWeight','bold','fontSize',20)
% legend('Trials', num2str(fit1(2)));
% ylabel('MPG','fontWeight','bold','fontSize',20)
% xlabel('Delta SOC','fontWeight','bold','fontSize',20);
% set(gca,'fontSize',20,'fontWeight','bold'),grid