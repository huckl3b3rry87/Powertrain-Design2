clear
clc
load small_parallel_HWFET_results

Results_Improvement = zeros(1,3);
Results_Pars = zeros(1,3);

figure(1);
mpg = No_Shift_UDDS_optimization.GLOBAL.f_min_hist.c(:,3);
iter = No_Shift_UDDS_optimization.GLOBAL.f_min_hist.iter;
plot(iter,mpg,'linewidth',8)
xlabel('Iterations')
ylabel('MPG'),grid
%title('No Shift UDDS')
set(gca,'FontSize',15,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',22,'fontWeight','bold')

figure(2);
CO = No_Shift_UDDS_optimization.GLOBAL.f_min_hist.c(:,4);
iter = No_Shift_UDDS_optimization.GLOBAL.f_min_hist.iter;
plot(iter,CO,'linewidth',8)
xlabel('Iterations')
ylabel('Carbon Monoxide (g)'),grid
%title('No Shift UDDS')
set(gca,'FontSize',15,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',22,'fontWeight','bold')

figure(3);
HC = No_Shift_UDDS_optimization.GLOBAL.f_min_hist.c(:,5);
iter = No_Shift_UDDS_optimization.GLOBAL.f_min_hist.iter;
plot(iter,HC,'linewidth',8)
xlabel('Iterations')
ylabel('Hydrocarbons (g)'),grid
%title('No Shift UDDS')
set(gca,'FontSize',15,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',22,'fontWeight','bold')

figure(4);
NOx = No_Shift_UDDS_optimization.GLOBAL.f_min_hist.c(:,6);
iter = No_Shift_UDDS_optimization.GLOBAL.f_min_hist.iter;
plot(iter,NOx,'linewidth',8)
xlabel('Iterations')
ylabel('Oxides of Nitrogen (g)'),grid
%title('No Shift UDDS')
set(gca,'FontSize',15,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',22,'fontWeight','bold')

figure(5);
PM = No_Shift_UDDS_optimization.GLOBAL.f_min_hist.c(:,7);
iter = No_Shift_UDDS_optimization.GLOBAL.f_min_hist.iter;
plot(iter,PM,'linewidth',8)
xlabel('Iterations')
ylabel('Particulate Matter (g)'),grid
%title('No Shift UDDS')
set(gca,'FontSize',15,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',22,'fontWeight','bold')

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

fprintf('\nPM\n')
First = PM(1)
Final = PM(length(PM))
Percent_improvement = (PM(length(PM))-PM(1))/PM(length(PM))*100
Results_Improvement(5,:) = [First, Final, Percent_improvement];

% Pull out Parameters
trq_eng = No_Shift_UDDS_optimization.GLOBAL.f_min_hist.x(:,1);
trq_mot = No_Shift_UDDS_optimization.GLOBAL.f_min_hist.x(:,2);
modules = No_Shift_UDDS_optimization.GLOBAL.f_min_hist.x(:,3);
fd_ratio = No_Shift_UDDS_optimization.GLOBAL.f_min_hist.x(:,4);
cap_scale = No_Shift_UDDS_optimization.GLOBAL.f_min_hist.x(:,5);
M_ratio = No_Shift_UDDS_optimization.GLOBAL.f_min_hist.x(:,6);
E_ratio = No_Shift_UDDS_optimization.GLOBAL.f_min_hist.x(:,7);

x_L=[0.05*trq_eng,0.05*trq_mot,0.05*modules,0.05*fd_ratio,0.05*cap_scale,0.05*M_ratio,0.05*E_ratio,0.05*Eng_weight];
x_U=[1.95*trq_eng,1.95*trq_mot,1.95*modules,1.95*fd_ratio,1.95*cap_scale,1.95*M_ratio,1.95*E_ratio,1.95*Eng_weight];

Results_Pars(1,:) = [x_L(1), x_L(1), trq_eng(end)];
Results_Pars(2,:) = [x_L(2), x_L(2), trq_mot(end)];
Results_Pars(3,:) = [x_L(3), x_L(3), modules(end)];
Results_Pars(4,:) = [x_L(4), x_L(4), fd_ratio(end)];
Results_Pars(5,:) = [x_L(5), x_L(5), cap_scale(end)];
Results_Pars(6,:) = [x_L(6), x_L(6), M_ratio(end)];
Results_Pars(7,:) = [x_L(7), x_L(7), E_ratio(end)];

if (size(No_Shift_UDDS_optimization.GLOBAL.f_min_hist.x,2) == 8)
    Eng_weight = No_Shift_UDDS_optimization.GLOBAL.f_min_hist.x(:,8);
    Results_Pars(8,:) = [x_L(8), x_L(8), Eng_weight(end)];
end

% Add in Grade Test, Accel Test Results, delta soc, and delta trace
time0_60 = No_Shift_UDDS_optimization.GLOBAL.f_min_hist.c(:,8);
time40_60 = No_Shift_UDDS_optimization.GLOBAL.f_min_hist.c(:,9);
time0_85 = No_Shift_UDDS_optimization.GLOBAL.f_min_hist.c(:,10);
grade = No_Shift_UDDS_optimization.GLOBAL.f_min_hist.c(:,11);
delta_soc = No_Shift_UDDS_optimization.GLOBAL.f_min_hist.c(:,1);
delta_trace = No_Shift_UDDS_optimization.GLOBAL.f_min_hist.c(:,2);

Results_Improvement(6,:) = [time0_60(1), time0_60(end), (time0_60(end)-time0_60(1))/time0_60(end)*100];
Results_Improvement(7,:) = [time40_60(1), time40_60(end), (time40_60(end)-time40_60(1))/time40_60(end)*100];
Results_Improvement(8,:) = [time0_85(1), time0_85(end), (time0_85(end)-time0_85(1))/time0_85(end)*100];
Results_Improvement(8,:) = [grade(1), grade(end), (grade(end)-grade(1))/grade(end)*100];
Results_Improvement(8,:) = [delta_soc(1), delta_soc(end), (delta_soc(end)-delta_soc(1))/delta_soc(end)*100];
Results_Improvement(8,:) = [delta_trace(1), delta_trace(end), (delta_trace(end)-delta_trace(1))/delta_trace(end)*100];

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

fprintf('\nBattery Modules')
First = ((8.3645^2)/0.0377)*modules(1)*cap_scale(1)
Final = ((8.3645^2)/0.0377)*modules(end)*cap_scale(end)
Increase = (Final-First)/Final*100
Results_Improvement(11,:) = [First, Final, Increase];

% Calculate the Mass of the Components

% Engine
First = (4.1384*trq_eng(1)*67000/1000+14.5);
Final = (4.1384*trq_eng(end)*67000/1000+14.5);
Increase = (Final-First)/Final*100;
Results_Improvement(12,:) = [First, Final, Increase];
% Motor
First = (0.8635*trq_mot(1)*49000/1000+28.5);
Final = (0.8635*trq_mot(end)*49000/1000+28.5);
Increase = (Final-First)/Final*100;
Results_Improvement(13,:) = [First, Final, Increase];
% Battery
First = ((44*.4536)/20)*cap_scale(1)*modules(1);
Final = ((44*.4536)/20)*cap_scale(end)*modules(end);
Increase = (Final-First)/Final*100;
Results_Improvement(14,:) = [First, Final, Increase];


