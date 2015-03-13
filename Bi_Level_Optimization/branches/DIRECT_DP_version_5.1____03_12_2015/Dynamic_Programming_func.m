% Dynamic Programming - written by Huckleberry Febbo - 07/20/2014

function [FAIL, MPG, emission, delta_SOC, sim] = Dynamic_Programming_func(param, vinf, dvar, cyc_data, RUN_TYPE, weight)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%---------------------Create New Tables-----------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

tables = [cyc_data.cyc_name, ' TABLES'];
mkdir(tables)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%---------------------Simulate all Possible Dynamics----------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%% Define State Grids
x1_grid = [0.4:RUN_TYPE.soc_size:0.8]';       % SOC
x1_length = length(x1_grid);

x2_grid = [0 1];   % Engine off (0)   &   Engine on (1)
x2_length = length(x2_grid);

x3_grid = vinf.gear;  % [1st 2nd...]             % Gear Level
x3_length = length(x3_grid);

% Define Control Grids
u1_grid = vinf.eng_control_trq';          % Engine Torque [N-m]
u1_length = length(u1_grid);

u2_grid = [-1 0 1];             % Shift down, nothing, up
u2_length = length(u2_grid);

%% Define Some soft cnstraints stuff
SOC_penalty = linspace(0.1,10,20);
NEAR_SOC_min = param.MIN_SOC + fliplr(linspace(0.001,0.02,20));
NEAR_SOC_max = param.MAX_SOC - fliplr(linspace(0.001,0.02,20));

% Make a New Folder
rmdir(tables,'s')                 % Delete any left over info
tables = [cyc_data.cyc_name, ' TABLES'];
mkdir(tables)
cd(tables);

if RUN_TYPE.sim == 0
    tic
end

for t = 1:cyc_data.time_cyc
    Tm_max = single(zeros(x2_length,x3_length,u1_length,u2_length));               % [x2]x[u1]x[u2]
    Tm_min = single(zeros(x2_length,x3_length,u1_length,u2_length));               % [x2]x[u1]x[u2]
    Tm_save = single(zeros(x2_length,x3_length,u1_length,u2_length));              % [x2]x[u1]x[u2]
    Wm_save = single(zeros(x2_length,x3_length,u1_length,u2_length));              % [x2]x[u1]x[u2]
    We_save = single(zeros(x2_length,x3_length,u1_length,u2_length));              % [x2]x[u1]x[u2]
    table_x1 = single(zeros(x1_length,x2_length,x3_length,u1_length,u2_length));   % [x2]x[u1]x[u2]
    inst_fuel = single(zeros(x2_length,x3_length,u1_length,u2_length));            % [x2]x[u1]x[u2]
    infeasible_Te = single(zeros(x2_length,x3_length,u1_length,u2_length));        % [x2]x[u1]x[u2]
    infeasible_Pbatt = single(zeros(x1_length,x2_length,x3_length,u1_length,u2_length));
    
    for x3 = 1:x3_length           % Go through all of the gears
        x3_c = x3_grid(x3);
        
        for u2 = 1:u2_length       % Shift down, don't shift, and shift up
            u2_c = u2_grid(u2);
            
            if (x3_c == x3_grid(1) && u2_c == u2_grid(1)) || (x3_c == x3_grid(x3_length) && u2_c == u2_grid(u2_length))
                u2_c = 0;          % Cannot shift
            end
            
            if u2 == 1 || u2 == 3  % Shift penalty
                Shift_Penalty = weight.shift;
            else
                Shift_Penalty = 0;
            end
            
            % Update x3 and Define New Gear ID
            New_Gear_Index = x3 + u2_c;
            x3_n = x3_grid(New_Gear_Index);
            
            for x2 = 1:x2_length             % Current Engine State
                ENG_state_c = x2_grid(x2);
                
                if  ENG_state_c == 0 
                    Eng_Penalty = [0; weight.engine_event*ones((length(u1_grid)-1),1)];  % [u1}x[1]
                else
                    Eng_Penalty = zeros(size(u1_grid));   % [u1}x[1]
                end
                We_c = cyc_data.Ww(t)*dvar.FD*x3_n;                % [rad/sec]
                Te_c =  u1_grid;                      % Engine Control, [u1]x[1]
                Te_drive = u1_grid;  % could be a difference form the auxillary power
                Wm_c = cyc_data.Ww(t)*dvar.FD*dvar.G;              % [rad/sec]
                Tm_c = cyc_data.Tw(t)/(dvar.FD*dvar.G)*ones(size(Te_c)) - Te_drive*x3_n/dvar.G;  % [u1]x[1]
                
                % Check Motor
                Tm_max_current = interp1(vinf.m_map_spd,vinf.m_max_trq,Wm_c)*ones(size(u1_grid));
                Tm_max(x2,x3,:,u2) =  Tm_max_current;    % [x2]x[x3]x[u1]x[u2]x[u3] 
                Tm_min_current = interp1(vinf.m_map_spd,vinf.m_max_gen_trq,Wm_c)*ones(size(u1_grid));
                Tm_min(x2,x3,:,u2) =  Tm_min_current;    % [x2]x[x3]x[u1]x[u2]x[u3]   
                
                Tm_save(x2,x3,:,u2) = Tm_c;                                                      % [x2]x[x3]x[u1]x[u2]x[u3]
                Wm_save(x2,x3,:,u2) = Wm_c*ones(size(u1_grid));                                  % [u1]x[1]  -  Check For Each Gear
                
                % Check Engine
                We_save(x2,x3,:,u2) = We_c*ones(size(u1_grid));                         % [u1]x[1]  - Check For Each Gear - Speed of the engine depends on the gear and engine state
                
                % Saturate Engine Speed - For Max Torque Lookup & Fuel Lookup
                if We_c < vinf.W_eng_min
                    We_c = vinf.W_eng_min;
                end
                if We_c > vinf.W_eng_max
                    We_c = vinf.W_eng_max;
                end
                
                Te_max =  interp1(vinf.eng_consum_spd_old,vinf.eng_max_trq,We_c)*ones(size(Te_c));
                infeasible_Te(x2,x3,:,u2) = (Te_max < Te_c);                            % [u1]x[1]
                
                Te_c(Te_c > Te_max) = Te_max(Te_c > Te_max);
                fuel = (interp2(vinf.eng_consum_trq',vinf.eng_consum_spd,vinf.eng_consum_fuel,Te_c,We_c,'linear')*cyc_data.dt)'; % sould modify maps near zero!!
                if RUN_TYPE.emiss == 1
                    NOx = (interp2(vinf.eng_consum_trq,vinf.eng_consum_spd,vinf.fc_nox_map,Te_c,We_c,'linear')*cyc_data.dt)';
                    CO = (interp2(vinf.eng_consum_trq,vinf.eng_consum_spd,vinf.fc_co_map,Te_c,We_c,'linear')*cyc_data.dt)';
                    HC = (interp2(vinf.eng_consum_trq,vinf.eng_consum_spd,vinf.fc_hc_map,Te_c,We_c,'linear')*cyc_data.dt)';
                end
                
                if RUN_TYPE.emiss == 1
                    inst_fuel(x2,x3,:,u2) = weight.fuel*fuel + weight.NOx*NOx + weight.CO*CO + weight.HC*HC + Shift_Penalty*ones(size(fuel))+ Eng_Penalty;
                else
                    inst_fuel(x2,x3,:,u2) = weight.fuel*fuel + Shift_Penalty*ones(size(fuel))+ Eng_Penalty;
                end
                
                % Update x1
                % Saturate the motor for the efficiency lookup table
                Tm_eff = Tm_c;
                Tm_eff(Tm_c > Tm_max_current) = Tm_max_current(Tm_c > Tm_max_current);
                Tm_eff(Tm_c < Tm_min_current) =  Tm_min_current(Tm_c < Tm_min_current);
                Wm_eff = Wm_c;
                Wm_eff(Wm_c > vinf.Wm_max) = vinf.Wm_max;
                Wm_eff(Wm_c < vinf.Wm_min) = vinf.Wm_min;
                
                eff_m = interp2(vinf.m_map_trq, vinf.m_map_spd, vinf.m_eff_map, Tm_eff, abs(Wm_eff))';
                eff_m(isnan(eff_m)) = 0.2;
                
                Pbat_charge = (Wm_c*Tm_c).*(eff_m*vinf.ess_coulombic_eff);    % Tm_c < 0
                Pbat_discharge = (Wm_c*Tm_c)./(eff_m*vinf.ess_coulombic_eff); % Battery needs to supply more power!!
                
                Pbat = Pbat_discharge;
                Pbat(Tm_c < 0) = Pbat_charge(Tm_c < 0);
                
                Pbat = repmat(Pbat,[1, x1_length]);
                Pbat = permute(Pbat, [2 1]);
                
                % Discharge
                Pbatt_max = repmat(interp1(vinf.ess_soc, vinf.ess_max_pwr_dis, x1_grid),[1,u1_length]);
                rint_discharge = repmat(interp1(vinf.ess_soc,vinf.ess_r_dis,x1_grid),[1,u1_length]);
                
                % Charge
                Pbatt_min = -repmat(interp1(vinf.ess_soc, vinf.ess_max_pwr_chg, x1_grid),[1,u1_length]);
                rint_charge = repmat(interp1(vinf.ess_soc,vinf.ess_r_chg,x1_grid),[1,u1_length]);
                
                % Check Battery Infeasibility
                infeasible_Pbatt(:,x2,x3,:,u2) = (Pbatt_max < Pbat); % Do not peanlize it for opperating too low, we can brake
                
                % Saturate the Battery
                Pbat(Pbatt_max < Pbat) = Pbatt_max(Pbatt_max < Pbat);
                Pbat(Pbat < Pbatt_min) = Pbatt_min(Pbat < Pbatt_min);
                
                % Charge & Discharge Resistances
                rint_c = rint_charge;
                rint_c(Pbat > 0) = rint_discharge(Pbat > 0);
                
                Voc_c = repmat(interp1(vinf.ess_soc,vinf.ess_voc,x1_grid), [1, u1_length]);
                SOC_c_matrix = repmat(x1_grid,[1, u1_length]);
                SOC_n =  SOC_c_matrix -(Voc_c -(Voc_c.^2 -4*Pbat.*rint_c).^(1/2))./(2*rint_c*vinf.ess_cap_ah*3600)*cyc_data.dt;
                
                table_x1(:,x2,x3,:,u2) = SOC_n;
                
            end       % End of x2 ( Engine State ) Loops
        end           % End of u2 ( Gear Control )
    end               % End of x3 (Gear State) loops
    % Check Motor
    infeasible_Tm = ((Tm_save > Tm_max)|(Tm_save < Tm_min));            % Can Brake to make the rest up
    infeasible_Tm = repmat(infeasible_Tm,[1,1,1,1,x1_length]);
    infeasible_Tm = permute(infeasible_Tm,[5 1 2 3 4]);
    
    infeasible_Wm = (Wm_save > vinf.Wm_max) | (Wm_save < vinf.Wm_min);              % [x2]x[x3]x[u1]x[u2]x[u3]
    infeasible_Wm = repmat(infeasible_Wm,[1,1,1,1,x1_length]);
    infeasible_Wm = permute(infeasible_Wm,[5 1 2 3 4]);
    
    % Check Engine
    infeasible_We = (We_save < 550*param.rpm2rads) | (We_save > vinf.W_eng_max);
%     infeasible_We(:,1,:,:)= zeros(size(infeasible_We(:,1,:,:)));  % Do not penalize the engine speed if ther eng is off
    infeasible_We = repmat(infeasible_We,[1,1,1,1,x1_length]);
    infeasible_We = permute(infeasible_We,[5 1 2 3 4]);
    infeasible_We(:,1,:,1,:)= zeros(size(infeasible_We(:,1,:,1,:)));  % Do not penalize the engine speed if ther eng is off
    
    infeasible_Te = repmat(infeasible_Te,[1,1,1,1,x1_length]);
    infeasible_Te = permute(infeasible_Te,[5 1 2 3 4]);
    
    % Check SOC
    infeasible_SOC = (table_x1 < param.MIN_SOC) | (table_x1 > param.MAX_SOC);        % [x2]x[x3]x[u1]x[u2]x[u3]
    table_x1(table_x1 > param.MAX_SOC) = param.MAX_SOC;
    table_x1(table_x1 < param.MIN_SOC) = param.MIN_SOC;
    
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
    SOC_soft = SOC_soft + SOC_penalty(20)*((param.MIN_SOC < table_x1) & (table_x1 < NEAR_SOC_min(20)));
    
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
    SOC_soft = SOC_soft + SOC_penalty(20)*((param.MAX_SOC > table_x1) & (table_x1 > NEAR_SOC_max(20)));
    
    inst_fuel = repmat(inst_fuel,[1,1,1,1,x1_length]); % Add an extra dimension for the fuel table
    inst_fuel = permute(inst_fuel,[5 1 2 3 4]);
    
    table_L = inst_fuel + SOC_soft + 10*infeasible_SOC + weight.infeasible*(infeasible_We + infeasible_Tm + infeasible_Wm + infeasible_Te + infeasible_Pbatt);   %[x2]x[u1]x[u2]x[u3]
    
    savename = ['Transitional Cost = ',num2str(t),' Table.mat'];
    save(savename,'table_x1','table_L');
    
    if RUN_TYPE.sim == 0  % For DP only runs
        complete = 100  - (cyc_data.time_cyc - t)/(cyc_data.time_cyc)*100;
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
BETA = weight.CS;
Desired_SOC = 0.55; % When you do the optimization - Extract solutions from the middle

folder = [cyc_data.cyc_name, ' TABLES'];
cd(folder);      % Go into the tables that were previously made

J_STAR = repmat(BETA*(x1_grid - Desired_SOC).^2, [1, x2_length,x3_length]); % terminal state penalty, [x1]x[x2]
J_STAR(J_STAR~=0) = J_STAR(J_STAR~=0) + weight.SOC_final;


for t = cyc_data.time_cyc:-1:1
    loadfile_name = ['Transitional Cost = ',num2str(t),' TABLE.mat'];
    load(loadfile_name);
    SOC_State_Penalty = single(zeros(x1_length,x2_length,x3_length,u1_length,u2_length));
    for x2 = 1:x2_length
        for x3 = 1:x3_length
            for u1 = 1:u1_length
                for u2 = 1:u2_length
                    % Next Gear
                    if u2 == 1 && x3 == 1 || u2 == 3 && x3 == x3_length
                        u2_c = 0; % Cannot Shift
                        Infeasible_Shift = weight.infeasible*single(ones(x1_length,1,1,1));
                    else
                        u2_c = u2_grid(u2);
                        Infeasible_Shift = single(zeros(x1_length,1,1,1));
                    end
                    x3_n = x3 + u2_c;
                    
                    % Next Engine State
                    if  u1 == 1   % Engine is off
                        x2_n = 1;
                    else
                        x2_n = 2;  % eng is on
                    end
                    F  = griddedInterpolant(x1_grid,J_STAR(:,x2_n,x3_n),'linear');  % Penalizing where they land! - Just penalizing the SOC
                    SOC_State_Penalty(:,x2,x3,u1,u2) = F(table_x1(:,x2,x3,u1,u2)) + Infeasible_Shift;
                end
            end
        end
    end
    J_temp = table_L + SOC_State_Penalty;
    
    for x1 = 1:x1_length
        for x2 = 1:x2_length
            for x3 = 1:x3_length
                S = squeeze(J_temp(x1,x2,x3,:,:));
                [minS,idx] = min(S(:));
                [u1_id,u2_id] = ind2sub(size(S),idx);
                opt_trq(x1,x2,x3) = u1_grid(u1_id);
                opt_id_u2(x1,x2,x3) = u2_id;
                
                % Define the new optimum value
                opt_value(x1,x2,x3) = J_temp(x1,x2,x3,u1_id,u2_id);  % Using Optimium Control Sequence [u1opt,u2opt,u3opt]
            end
        end
    end
    J_STAR = opt_value;   % Next terminal cost!
    
    savename=['Cost & Control = ',num2str(t),' TABLE.mat'];
    save(savename,'J_STAR','opt_trq','opt_id_u2');
    
    if RUN_TYPE.sim == 0
        complete = (cyc_data.time_cyc - t)/(cyc_data.time_cyc)*100;
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
%--------------------------Simulate Final Run-----------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
if RUN_TYPE.sim == 0
    clc
    min_t  = toc/60;
    fprintf('__________________________________________________________________\n\n')
    fprintf('Total Time (min) for Dynamics and Dynamic Programming = ')
    fprintf(num2str(min_t))
    fprintf('\n')
    fprintf('__________________________________________________________________\n\n')
end

% Pre-Allocate for Speed
fail_inner_SOC = zeros(1,cyc_data.time_cyc);
fail_inner_Te = zeros(1,cyc_data.time_cyc);
fail_inner_We = zeros(1,cyc_data.time_cyc);
fail_inner_Tm = zeros(1,cyc_data.time_cyc);
fail_inner_Wm = zeros(1,cyc_data.time_cyc);
fail_inner_Pbatt = zeros(1,cyc_data.time_cyc);
fail_inner_Shift = zeros(1,cyc_data.time_cyc);

if RUN_TYPE.sim == 0
    fprintf('-------------------------------------------------\n');
    fprintf('Simulating Final Run! ')
    fprintf('\n')
    fprintf('-------------------------------------------------\n');
end

folder = [cyc_data.cyc_name, ' TABLES'];
cd(folder);

% Pre-Allocate for Speed
sim.SOC_final = zeros(1,cyc_data.time_cyc);
sim.inst_fuel = zeros(1,cyc_data.time_cyc);
sim.W_mot = zeros(1,cyc_data.time_cyc);
sim.T_mot = zeros(1,cyc_data.time_cyc);
sim.W_eng = zeros(1,cyc_data.time_cyc);
sim.T_eng = zeros(1,cyc_data.time_cyc);
sim.ENG_state = zeros(1,cyc_data.time_cyc);
sim.GEAR = zeros(1,cyc_data.time_cyc);
sim.U1_sim = zeros(1,cyc_data.time_cyc);
sim.U2_sim = zeros(1,cyc_data.time_cyc);
sim.GEAR_save = zeros(1,cyc_data.time_cyc);
sim.J = zeros(1,cyc_data.time_cyc);
sim.Pbatt_sim = zeros(1,cyc_data.time_cyc);

% Define Initial Conditions
SOC_c = 0.55;
x2 = 1;                    % Engine Off
x3 = 1;                    % Start in First Gear

for t = 1:1:cyc_data.time_cyc
   % Save States For Simulation Results
    sim.SOC_final(t) = SOC_c;
    sim.ENG_state(t) = x2-1;    % Just on or off  (0 or 1)
    sim.GEAR(t) = x3;
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    %------------- Load & Determine all Control Signals --------------%
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    load(['Cost & Control = ',num2str(t),' TABLE.mat']);
    
    % Engine Torque Control
    Te_c = interp1(x1_grid,opt_trq(:,x2,x3),SOC_c,'linear');
    if x2 == 1 && Te_c ~= 0; ENG_event = 1; else ENG_event = 0; end
    
    % Shifting Control
    id_lookup_u2 = interp1(x1_grid,opt_id_u2(:,x2,x3),SOC_c,'nearest');  %Use index extraction!!
    u2_c = u2_grid(id_lookup_u2);
    if u2_c == 0; SHIFT_event = 0; else SHIFT_event = 1; end
    
    sim.J(t) =  interp1(x1_grid,J_STAR(:,x2,x3),SOC_c,'linear'); 
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    
    % Shifting Control Check
    if (x3 == 1 && u2_c == u2_grid(1)) || (x3 == x3_length && u2_c == u2_grid(u2_length))
        FAIL_Shift = 1;
    else
        FAIL_Shift = 0;
    end
    
    % Update x3 and Define New Gear ID
    New_Gear_Index = x3 + u2_c;
    x3_n = x3_grid(New_Gear_Index);
    %
    We_c = cyc_data.Ww(t)*dvar.FD*x3_n;                % [rad/sec]
    Te_drive = Te_c;
    Wm_c = cyc_data.Ww(t)*dvar.FD*dvar.G;              % [rad/sec]
    Tm_c = cyc_data.Tw(t)/(dvar.FD*dvar.G)*ones(size(Te_c)) - Te_drive*x3_n/dvar.G;  % [u1]x[1]
    
    if We_c < vinf.W_eng_min
        We_fuel = vinf.W_eng_min;
    elseif We_c > vinf.W_eng_max
        We_fuel = vinf.W_eng_max;
    else

        We_fuel = We_c;
    end
    
    if (Te_c ~= 0 && ((We_c < 550*param.rpm2rads) || We_c > vinf.W_eng_max))
        Fail_We = 1;
    else
        Fail_We = 0;
    end
    %
    Te_max =  interp1(vinf.eng_consum_spd_old,vinf.eng_max_trq,We_fuel); % Canged tp the saturated value
    
    if Te_c > Te_max;
        Fail_Te = 1;
        Te_fuel = Te_max;
    else
        Fail_Te = 0;
        Te_fuel = Te_c;
    end
    fuel = (interp2(vinf.eng_consum_trq',vinf.eng_consum_spd,vinf.eng_consum_fuel,Te_fuel,We_fuel,'linear')*cyc_data.dt)'; % sould modify maps near zero!!
    if RUN_TYPE.emiss == 1
        NOx = (interp2(vinf.eng_consum_trq,vinf.eng_consum_spd,vinf.fc_nox_map,Te_fuel,We_fuel,'linear')*cyc_data.dt)';
        CO = (interp2(vinf.eng_consum_trq,vinf.eng_consum_spd,vinf.fc_co_map,Te_fuel,We_fuel,'linear')*cyc_data.dt)';
        HC = (interp2(vinf.eng_consum_trq,vinf.eng_consum_spd,vinf.fc_hc_map,Te_fuel,We_fuel,'linear')*cyc_data.dt)';
    end
    
    %             % Update x1
    % Saturate the motor for the efficiency lookup table
    Tm_max= interp1(vinf.m_map_spd,vinf.m_max_trq,Wm_c);
    if Tm_c > Tm_max;
        Fail_Tm = 1;
        Tm_eff = Tm_max;
    elseif Tm_c < -Tm_max;
        Fail_Tm = 1;
        Tm_eff = -Tm_max;
    else
        Fail_Tm = 0;
        Tm_eff = Tm_c;
    end
    
    if Wm_c > vinf.Wm_max
        Fail_Wm = 1;
        Wm_eff = vinf.Wm_max;
    elseif Wm_c < vinf.Wm_min
        Fail_Wm = 1;
        Wm_eff = vinf.Wm_min;
    else
        Fail_Wm = 0;
        Wm_eff = Wm_c;
    end
    
    eff_m = interp2(vinf.m_map_trq, vinf.m_map_spd, vinf.m_eff_map, Tm_eff, abs(Wm_eff))';
    
    if isnan(eff_m)
        eff_m = 0.2;
    end
    
    Pbat_charge = (Wm_c*Tm_c).*(eff_m*vinf.ess_coulombic_eff);    % Tm_c < 0            
    Pbat_discharge = (Wm_c*Tm_c)./(eff_m*vinf.ess_coulombic_eff); % Battery needs to supply more power!!
    
    %     Pbat = Pbat_discharge;
    %     Pbat(Tm_c < 0) = Pbat_charge(Tm_c < 0);  % again should not be saturated one here!!
    if Tm_c < 0;
        Pbat = Pbat_charge;
    else
        Pbat = Pbat_discharge;
    end
    
    
    % Discharge
    Pbatt_max = interp1(vinf.ess_soc, vinf.ess_max_pwr_dis, SOC_c);
    rint_discharge = interp1(vinf.ess_soc,vinf.ess_r_dis,SOC_c);
    
    % Charge
    Pbatt_min = -interp1(vinf.ess_soc, vinf.ess_max_pwr_chg,SOC_c);
    rint_charge = interp1(vinf.ess_soc,vinf.ess_r_chg,SOC_c);
    
    % Saturate the Battery
    if Pbat > Pbatt_max
        FAIL_Pbatt = 1;
        Pbat_eff = Pbatt_max;
    elseif Pbat < Pbatt_min
        FAIL_Pbatt = 1;
        Pbat_eff = Pbatt_min;
    else
        FAIL_Pbatt = 0;
        Pbat_eff = Pbat;
    end
    
    % Charge & Discharge Resistances
    if Pbat > 0
        rint_c = rint_discharge;
    else
        rint_c = rint_charge;
    end
    
    Voc_c = interp1(vinf.ess_soc,vinf.ess_voc,SOC_c);
    SOC_n =  SOC_c -(Voc_c -(Voc_c.^2 -4*Pbat_eff.*rint_c).^(1/2))./(2*rint_c*vinf.ess_cap_ah*3600)*cyc_data.dt;
    
    % Check new SOC
    if SOC_n > param.MAX_SOC
        SOC_n = param.MAX_SOC;
        Fail_SOC = 1;
    elseif SOC_n < param.MIN_SOC
        SOC_n = param.MIN_SOC;
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
    fail_inner_Shift(t) = FAIL_Shift;
    
    % Save Remaining Simulation Variables
    if Te_c == 0; sim.W_eng(t) = 0; else sim.W_eng(t) = We_c; end
    sim.W_mot(t) = Wm_c;
    sim.T_eng(t) = Te_c;
    sim.T_mot(t) = Tm_c;
    sim.inst_fuel(t) = fuel;
    sim.SHIFT_Event(t) = SHIFT_event;
    sim.ENG_Event(t) = ENG_event;
    
    if RUN_TYPE.emiss == 1
        sim.NOx(t) = NOx;
        sim.CO(t) = CO;
        sim.HC(t) = HC;
    end
    sim.Pbatt_sim(t) = Pbat_eff;
    sim.eff_m_sim(t) = eff_m;
    
    % Plot the control signal indicies
    sim.U1_sim(t)  =  Tm_c;
    sim.U2_sim(t)  = id_lookup_u2;
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    %---------------------- Update the States-------------------------%
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    % Update x1 for next time step
    SOC_c = SOC_n;
    x3 = find(x3_n == x3_grid);
    % Update engine state for simulation
    if  Te_c == 0; x2 = 1; else x2 = 2; end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    
end

sim.T_eng(sim.T_eng<0) = 0;
delta_SOC = sim.SOC_final(end) - sim.SOC_final(1);
total_fuel_gram = sum(sim.inst_fuel); % dt is 1 - grams
sim.EE = sum(sim.ENG_Event);
sim.SE = sum(sim.SHIFT_Event);

if RUN_TYPE.emiss == 1
    emission.NOx = sum(sim.NOx);   % dt is 1 - grams
    emission.HC = sum(sim.HC);
    emission.CO = sum(sim.CO);
else
    emission = NaN;
end
total_distance_mile = sum(cyc_data.cyc_spd)/3600;
MPG = total_distance_mile/(total_fuel_gram/1000/param.gasoline_density*param.liter2gallon);
cd ..    % Come out of folder

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%-----------------------Check Outer Feasibility---------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%SOC
if any(any(fail_inner_SOC))
    FAIL.fail_outer_SOC = 1;
else
    FAIL.fail_outer_SOC = 0;
end
if abs(delta_SOC) > RUN_TYPE.soc_size
    FAIL.fail_dSOC = 1;
else
    FAIL.fail_dSOC = 0;
end
%We
if any(any(fail_inner_We))
    FAIL.fail_outer_We = 1;
else
    FAIL.fail_outer_We = 0;
end
%Te
if any(any(fail_inner_Te))
    FAIL.fail_outer_Te = 1;
else
    FAIL.fail_outer_Te = 0;
end
%Wm
if any(any(fail_inner_Wm))
    FAIL.fail_outer_Wm = 1;
else
    FAIL.fail_outer_Wm = 0;
end
%Tm
if any(any(fail_inner_Tm))
    FAIL.fail_outer_Tm = 1;
else
    FAIL.fail_outer_Tm = 0;
end
%Pbatt
if any(any(fail_inner_Pbatt))
    FAIL.fail_outer_Pbatt = 1;
else
    FAIL.fail_outer_Pbatt = 0;
end
%Shift
if any(any(fail_inner_Shift))
    FAIL.fail_outer_Shift = 1;
else
    FAIL.fail_outer_Shift = 0;
end
FAIL.final = ((FAIL.fail_outer_Tm + FAIL.fail_outer_Wm + FAIL.fail_outer_Te + FAIL.fail_outer_We + FAIL.fail_outer_SOC + FAIL.fail_dSOC + FAIL.fail_outer_Pbatt + FAIL.fail_outer_Shift) ~= 0);

end
