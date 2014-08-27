% Compare Drive Cycles
close all
load CYC_US06_AV;

AV_SPD = cyc_mph(:,2);

load CYC_US06;

SPD = cyc_mph(:,2);

plot(cyc_mph(:,1),SPD)
hold on
plot(cyc_mph(:,1),AV_SPD,'g')
legend('US06','Modified US06')
xlabel('time (s)')
ylabel('Speed (MPH)'),grid
% % Specify the position and the size of the 2. axis
% x_a = 225; y_a = 60; w_a = 30; h_a = 20;
% ax = axes('Units', 'Normalized', ...
%           'Position', [x_a, y_a, w_a, h_a], ...
%           'XTick', [], ...
%           'YTick', [], ...
%           'Box', 'on', ...
%           'LineWidth', 2, ...
%           'Color', [0.95, 0.99, 0.95]);
% hold on;
% % plot(x, y);
% xlabel('Detail at X==0.95');
% % axis([x_r-w_r/2, x_r+w_r/2, y_r-h_r/2, y_r+h_r/2]);
hold off

%% other cycles
clear
clc

load AA_final

cyc_mph(:,2) = smooth(cyc_mph(:,2), 'moving');

save('AA_final_AV','cyc_mph');

%% D10
clear
cyc_name = 'Natralistic Ann Arbor Drive Cycle';
load AA_final_AV;
AV_SPD = cyc_mph(:,2);

load AA_final;

SPD = cyc_mph(:,2);
figure(1);
plot(cyc_mph(:,1),SPD)
hold on
plot(cyc_mph(:,1),AV_SPD,'g')
legend('Actual Driver Data','Modified Driver Data')
ylabel('Speed (mph)','fontWeight','bold','fontSize',12)
xlabel('time (sec)','fontWeight','bold','fontSize',12);
set(gca,'fontSize',12,'fontWeight','bold'),grid on
title({cyc_name},'fontWeight','bold','fontSize',16)
xlim([0 length(cyc_mph(:,1))*1.5])
grid on;
magnifyOnFigure;

%% HWFET
clear
cyc_name = 'HWFET';
load CYC_HWFET_AV;

AV_SPD = cyc_mph(:,2);

load CYC_HWFET;

SPD = cyc_mph(:,2);
figure(2);
plot(cyc_mph(:,1),SPD)
hold on
plot(cyc_mph(:,1),AV_SPD,'g')
legend('HWFET','Modified HWFET')
ylabel('Speed (mph)','fontWeight','bold','fontSize',12)
xlabel('time (sec)','fontWeight','bold','fontSize',12);
set(gca,'fontSize',12,'fontWeight','bold')
title({cyc_name},'fontWeight','bold','fontSize',16)
grid on;
magnifyOnFigure;
