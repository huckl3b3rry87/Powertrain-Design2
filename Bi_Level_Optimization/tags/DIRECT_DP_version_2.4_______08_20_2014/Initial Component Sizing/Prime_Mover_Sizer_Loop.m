%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%------------------------Component Sizer----------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
clear all
close all
clc
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%---Determine Acceleration Req. Based off Drive Cycles--------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
S = 10;                    % Number of points to calcualte acceleration at
dt_2 = 0.1;
graph = 0;
[ V_0, V_f, Acc_Final ] = Accel_Req( S, dt_2, graph ); % V_Final is in MPH

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%------------------Load the Component Data--------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
cd('Components');

% Engine
Engine_41_kW;
% Engine_102_kW;
% Engine_2rz_0410;

% Motor
Motor_75_kW;

Battery_ADVISOR;

% Vehicle
Vehicle_Parameters_4_HI_AV;
cd ..

data;                              %  Put all the data into structures
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%-------------------------Design Variables--------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
dvar.fc_trq_scale = 1;  % Might need to run a loop on this guy - OR ~ Just play with it!!
dvar.mc_trq_scale = 1;
dvar.FD = 3;
dvar.G = 1;  % Also change this so that you can get to the maximum speed!
vinf.mc_max_pwr_kW =  dvar.mc_trq_scale*vinf.mc_max_pwr_kW;
dvar.module_number = ceil(4*mc_max_pwr_kW*1000*Rint_size/(Voc_size^2));
Manipulate_Data_Structure;     % Update the data based off of the new variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%--------------------------Engine Power-----------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
initial_fc_trq_scale = dvar.fc_trq_scale;
initial_FD = dvar.FD;

for X = 1:10
   
    test = 1; % Zero Grade (Max Speed)
    ttt = 1;
    for FD = 1:0.05:8
        dvar.FD = FD;
        [~,~, V_max_t, ~, V_max_6,~] = Engine_Power_Sizer( param, vinf, dvar, test );
        
        if isempty(V_max_6);
            V_max_6 = NaN;
        end
        
        if isempty(V_max_t)
            V_max_t = NaN;
        end
        
        FD_sim(ttt) = FD;
        V_6(ttt) = V_max_6;
        V_act(ttt) =  V_max_t;
        ttt = ttt + 1;
    end
    
    [Max_Spd, I] = max(V_6/param.mph_mps);
    Final_Drive_Plot;
    
    % Resimualte with selected gear ratio
    fprintf(' The Best Final Drive Ratio is: ')
    dvar.FD = FD_sim(I)
    [Sim_Grade,F_max_t, V_max_t, F_max_6, V_max_6, RR ] = Engine_Power_Sizer( param, vinf, dvar, test );
    Grade_Plot;
    
    % Simulation Results
     fc_trq(X) = dvar.fc_trq_scale;
     FD_save(X) = dvar.FD;
     clear I
     
    %% Now size the engine so that it meets the rest of the grade requirements
    test = 2;
    V_test = 55*param.mph_mps;
    iii = 1;
    for fc_trq_scale = 0.5:0.05:1.75
        dvar.fc_trq_scale = fc_trq_scale;
        Manipulate_Data_Structure;  % Update the data based off of the new motor power and battery module
        [Sim_Grade,F_max_t, V_max_t, F_max_6, V_max_6, RR ] = Engine_Power_Sizer( param, vinf, dvar, test );
        
        if isempty(V_max_t)
            V_diff(iii) = NaN;
        else
            V_diff(iii) = abs(V_max_t - V_test); % Needs to be greater!! -Change this later
        end
        fc_trq_sim(iii) =  fc_trq_scale;
        iii = iii + 1;
    end
    figure(3); clf
    plot(fc_trq_sim,V_diff/param.mph_mps)
    [Min, I] = min(V_diff);
    fprintf('The best fc_trq_scale is: ')
    fc_trq_sim(I)
    
    % Resimulate with selected trq scale
    dvar.fc_trq_scale = fc_trq_sim(I)

    Manipulate_Data_Structure;  % Update the data based off of the new motor power and battery module
    [Sim_Grade,F_max_t, V_max_t, F_max_6, V_max_6, RR ] = Engine_Power_Sizer( param, vinf, dvar, test );
    Grade_Plot;
   clear I
end
figure(3);
[AX,L1,L2] = plotyy(1:X, FD_save, 1:X,fc_trq);
set(L1,'marker','x', 'markersize', 5, 'markerf','b','linewidth',3);
set(L2,'marker','x', 'markersize', 5, 'markerf','g','linewidth',3);
set(AX, 'fontsize', 20,'fontweight','bold','XTickLabel',[]);
ylabel({'Final Drive Ratio'}, 'parent', AX(1));
ylabel('fc trq scale', 'parent', AX(2));grid
xlabel('Itterations')

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%---------------------------Motor Power-----------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%



% % After acceleration tests and dv have been updated
% Manipulate_Data_Structure;  % Update the data based off of the new motor power and battery module
%




