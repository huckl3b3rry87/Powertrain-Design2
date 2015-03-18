clear
close all
clc
load HWFET_50_1_Nm_005_SOC_no_emiss

Results_Improvement = zeros(1,3);
Results_Pars = zeros(1,3);

figure(1);
mpg = DP_optimization.GLOBAL.f_min_hist.c(:,3);
iter = DP_optimization.GLOBAL.f_min_hist.iter;
plot(iter,mpg,'linewidth',8)
xlabel('Iterations')
ylabel('MPG'),grid
%title('No Shift UDDS')
set(gca,'FontSize',15,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',22,'fontWeight','bold')



% Calculate the Improvement
fprintf('\nMPG\n')
First = mpg(1)
Final = mpg(length(mpg))
Percent_improvement = (mpg(length(mpg))-mpg(1))/mpg(length(mpg))*100
Results_Improvement(1,:) = [First, Final, Percent_improvement];


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


