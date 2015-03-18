clear;
clf;
close all;
C = jet;
h = 14;
r=6;
t = 10;
q = 20;  % For colors
start = 1;
stop = 10;
for i = [start, 2];
    if i == 1
        load HWFET_50_1_Nm_005_SOC_no_emiss;
        c = 'No Emissions';
    elseif i == 2
        load HWFET_50_1_Nm_005_SOC_emiss;
        c = 'Emissions';
    elseif i == 3
        load SOC_grid_005_trq_20_Nm_HWFET;
        c = '20 Nm';
    elseif i == 4
        load SOC_grid_005_trq_5_Nm_HWFET;
        c = '0.001';
    elseif i == 5
        load SOC_grid_0005;
        c = '0.0005';
    elseif i == 6
        load SOC_grid_005_interp;
        c = '0.005 lin. interp';
    elseif i == 7
        load SOC_grid_001_interp;
        c = '0.001 lin. interp.';
    elseif i == 8
        load SOC_grid_0005_interp;
        c = '0.0005 lin. interp.';
    elseif i == 9
        load SOC_grid_005_interp_cubic;
        c = '0.005 cubic interp.';
    elseif i == 10
        load SOC_grid_0005_interp_cubic;
        c = '0.0005 cubic. interp.';
    elseif i == 11;
        load SOC_grid_005_interp_full_cyc;
        c = '0.005 line. interp FULL';
    elseif i == 12;
        load SOC_grid_005_interp_lower_weights;
        c = '0.005 line. interp Lower weights';
    elseif i == 13;
        load SOC_grid_005_interp_full_cyc_smooth;
        c = '0.005 line. interp Full CYC Smoothed Data';
    else
        load SOC_grid_001_interp_full_cyc_smooth;
        c = '0.001 line. interp Full CYC Smoothed Data';
    end
    
    if i == start
        leg = {c};
        a = 1;
    else
        leg = {leg{1:a}, c};
        a = a+1;
    end
    
    iter = DP_optimization.GLOBAL.f_min_hist.iter;
    
    figure(2)
    mpg = DP_optimization.GLOBAL.f_min_hist.c(:,3);
    plot(iter,mpg,'-ko',...
        'LineWidth',0.5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',C(a+(a-1)*q,:),...
        'MarkerSize',r)
    grid on
    ylabel('MPG')
    xlabel('Itterations')
    legend(leg{1:a})
    set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','bold')
    hold on
    
    figure(3)
    CO = DP_optimization.GLOBAL.f_min_hist.c(:,5);
    plot(iter,CO,'-ko',...
        'LineWidth',0.5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',C(a+(a-1)*q,:),...
        'MarkerSize',r)
    grid on
    ylabel('CO (g)')
    xlabel('Itterations')
    legend(leg{1:a})
    set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','bold')
    hold on
    
    figure(4)
    HC = DP_optimization.GLOBAL.f_min_hist.c(:,6);
    plot(iter,HC,'-ko',...
        'LineWidth',0.5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',C(a+(a-1)*q,:),...
        'MarkerSize',r)
    grid on
    ylabel('HC (g)')
    xlabel('Itterations')
    legend(leg{1:a})
    set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','bold')
    hold on
    
    figure(5)
    NOx = DP_optimization.GLOBAL.f_min_hist.c(:,4);
    plot(iter,NOx,'-ko',...
        'LineWidth',0.5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',C(a+(a-1)*q,:),...
        'MarkerSize',r)
    grid on
    ylabel('NOx')
    xlabel('Iterations')
    legend(leg{1:a})
    set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','bold')
    hold on
    
end

hold off
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','bold')



Results_Improvement = zeros(1,3);
Results_Pars = zeros(1,3);
% Calculate the Improvement
fprintf('\nMPG\n')
First = mpg(1)
Final = mpg(length(mpg))
Percent_improvement = (mpg(length(mpg))-mpg(1))/mpg(length(mpg))*100
Results_Improvement(1,:) = [First, Final, Percent_improvement];

fprintf('\nCO\n')
First = CO(1)
Final = CO(length(CO))
Percent_improvement = (CO(length(CO))-CO(1))/CO(length(CO))*100
Results_Improvement(2,:) = [First, Final, Percent_improvement];

fprintf('\nHC\n')
First = HC(1)
Final = HC(length(HC))
Percent_improvement = (HC(length(HC))-HC(1))/HC(length(HC))*100
Results_Improvement(3,:) = [First, Final, Percent_improvement];

fprintf('\nNOx\n')
First = NOx(1)
Final = NOx(length(NOx))
Percent_improvement = (NOx(length(NOx))-NOx(1))/NOx(length(NOx))*100
Results_Improvement(4,:) = [First, Final, Percent_improvement];


% Pull out Parameters
trq_eng = DP_optimization.GLOBAL.f_min_hist.x(:,3);
trq_mot = DP_optimization.GLOBAL.f_min_hist.x(:,4);
fd_ratio = DP_optimization.GLOBAL.f_min_hist.x(:,1);
G_ratio = DP_optimization.GLOBAL.f_min_hist.x(:,2);


x_L=[0.5*trq_eng,0.5*trq_mot,0.5*fd_ratio,0.5*G_ratio];
x_U=[1.5*trq_eng,1.5*trq_mot,1.5*fd_ratio,1.5*G_ratio];

Results_Pars(1,:) = [x_L(1), x_L(1), trq_eng(end)];
Results_Pars(2,:) = [x_L(2), x_L(2), trq_mot(end)];
Results_Pars(3,:) = [x_L(3), x_L(3), fd_ratio(end)];
Results_Pars(4,:) = [x_L(4), x_L(4), G_ratio(end)];


% % Add in Grade Test, Accel Test Results, delta soc, and delta trace
% time0_60 = DP_optimization.GLOBAL.f_min_hist.c(:,8);
% time40_60 = DP_optimization.GLOBAL.f_min_hist.c(:,9);
% time0_85 = DP_optimization.GLOBAL.f_min_hist.c(:,10);
% grade = DP_optimization.GLOBAL.f_min_hist.c(:,11);
% delta_soc = DP_optimization.GLOBAL.f_min_hist.c(:,1);
% delta_trace = DP_optimization.GLOBAL.f_min_hist.c(:,2);
% 
% Results_Improvement(6,:) = [time0_60(1), time0_60(end), (time0_60(end)-time0_60(1))/time0_60(end)*100];
% Results_Improvement(7,:) = [time40_60(1), time40_60(end), (time40_60(end)-time40_60(1))/time40_60(end)*100];
% Results_Improvement(8,:) = [time0_85(1), time0_85(end), (time0_85(end)-time0_85(1))/time0_85(end)*100];
% Results_Improvement(8,:) = [grade(1), grade(end), (grade(end)-grade(1))/grade(end)*100];
% Results_Improvement(8,:) = [delta_soc(1), delta_soc(end), (delta_soc(end)-delta_soc(1))/delta_soc(end)*100];
% Results_Improvement(8,:) = [delta_trace(1), delta_trace(end), (delta_trace(end)-delta_trace(1))/delta_trace(end)*100];

% Calcualte the Size of the Prime Movers
fprintf('\nPrime Mover Sizes\n')

fprintf('\nICE')
First = 67*trq_eng(1)
Final = 67*trq_eng(length(trq_eng))
Increase = (Final-First)/Final*100
Results_Improvement(9,:) = [First, Final, Increase];

fprintf('\nMotor')
First = 49*trq_mot(1)
Final = 49*trq_mot(length(trq_mot))
Increase = (Final-First)/Final*100
Results_Improvement(10,:) = [First, Final, Increase];



% % Calculate the Mass of the Components
% 
% % Engine
% First = (4.1384*trq_eng(1)*67000/1000+14.5);
% Final = (4.1384*trq_eng(end)*67000/1000+14.5);
% Increase = (Final-First)/Final*100;
% Results_Improvement(12,:) = [First, Final, Increase];
% % Motor
% First = (0.8635*trq_mot(1)*49000/1000+28.5);
% Final = (0.8635*trq_mot(end)*49000/1000+28.5);
% Increase = (Final-First)/Final*100;
% Results_Improvement(13,:) = [First, Final, Increase];
% % Battery
% First = ((44*.4536)/20)*cap_scale(1)*modules(1);
% Final = ((44*.4536)/20)*cap_scale(end)*modules(end);
% Increase = (Final-First)/Final*100;
% Results_Improvement(14,:) = [First, Final, Increase];



