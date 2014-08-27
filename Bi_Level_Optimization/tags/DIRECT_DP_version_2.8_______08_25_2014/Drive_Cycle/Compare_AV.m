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

load D10

cyc_mph(:,2) = smooth(cyc_mph(:,2), 'moving');

save('D10AV','cyc_mph');

%% D10
clear
cyc_name = 'Natralistic Ann Arbor Drive Cycle';
load D10AV;
AV_SPD = cyc_mph(:,2)*2.23694;

load D10;

SPD = cyc_mph(:,2)*2.23694;
figure(1);
plot(cyc_mph(:,1),SPD)
hold on
plot(cyc_mph(:,1),AV_SPD,'g')
legend('Driver 10','Modified Driver 10')
ylabel('Speed (mph)','fontWeight','bold','fontSize',12)
xlabel('time (sec)','fontWeight','bold','fontSize',12);
set(gca,'fontSize',12,'fontWeight','bold'),grid on
title({cyc_name},'fontWeight','bold','fontSize',16)
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
%%
clear
cyc_name = 'Natralistic Ann Arbor Drive Cycle';
dt = 1;
load D8;
SPD = cyc_mph(:,2)*2.23694; % Assuming data is in m/s convert to mph 
v = cyc_mph(:,2);
% SPD = cyc_mph(:,2);
figure(1);
plot(cyc_mph(:,1),SPD)

legend('Driver 10')
ylabel('Speed (mph)','fontWeight','bold','fontSize',12)
xlabel('time (sec)','fontWeight','bold','fontSize',12);
set(gca,'fontSize',12,'fontWeight','bold'),grid on
title({cyc_name},'fontWeight','bold','fontSize',16)
grid on;
for i = 2:size(SPD)   % Backwards difference omits first point
    a(i) = (v(i) - v(i-1))/dt;
end
figure(2)
plot(cyc_mph(:,1),a)
legend('Acceleration (m/s^2)')

%%
clear
mph_mps = 0.44704; 
cyc_name = 'Natralistic Ann Arbor Drive Cycle';
dt = 1;
load D6;
SPD = cyc_mph(:,2); 
v = cyc_mph(:,2)*mph_mps; % Assuming data is in mph convert to m/s
% SPD = cyc_mph(:,2);
figure(1);
plot(cyc_mph(:,1),SPD)

legend('Driver 10')
ylabel('Speed (mph)','fontWeight','bold','fontSize',12)
xlabel('time (sec)','fontWeight','bold','fontSize',12);
set(gca,'fontSize',12,'fontWeight','bold'),grid on
title({cyc_name},'fontWeight','bold','fontSize',16)
grid on;
for i = 2:size(SPD)   % Backwards difference omits first point
    a(i) = (v(i) - v(i-1))/dt;
end
figure(2)
plot(cyc_mph(:,1),a)
legend('Acceleration (m/s^2)')

%%
clear
kmh_mph = 0.621; 
cyc_name = 'Natralistic Ann Arbor Drive Cycle';
dt = 1;
load D2;
SPD = cyc_mph(:,2)*kmh_mph; 
v = cyc_mph(:,2)*0.27777; % Assuming data is in mph convert to m/s
% SPD = cyc_mph(:,2);
figure(1);
plot(cyc_mph(:,1),SPD)

legend('Driver 10')
ylabel('Speed (mph)','fontWeight','bold','fontSize',12)
xlabel('time (sec)','fontWeight','bold','fontSize',12);
set(gca,'fontSize',12,'fontWeight','bold'),grid on
title({cyc_name},'fontWeight','bold','fontSize',16)
grid on;
for i = 2:size(SPD)   % Backwards difference omits first point
    a(i) = (v(i) - v(i-1))/dt;
end
figure(2)
plot(cyc_mph(:,1),a)
legend('Acceleration (m/s^2)')