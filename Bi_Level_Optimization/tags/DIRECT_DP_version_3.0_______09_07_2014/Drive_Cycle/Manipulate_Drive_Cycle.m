% cyc_name = 'HWFET';
% load CYC_HWFET
cyc_name = 'AA Final';
load AA_final

% Modify The Drive Cycle For DP
dt = 1;

v = cyc_mph(:,2)*mph_mps;  % Cycle Speed in (m/s)
time_cyc = length(v);
cyc_time = 1:time_cyc;
% Need to define the a(t) of the profile
for i = 2:size(v)   % Backwards difference omits first point
    a(i) = (v(i) - v(i-1))/dt;
end

for t = 1:1:time_cyc  

    %% Calculate Tw, Ww and Pd from the Drive Cycle
    if v(t) ~= 0
        Tw(t) = vinf.rwh*(vinf.m*param.g*sin(param.grade) + vinf.Frr*vinf.m*param.g*cos(param.grade) + 0.5*param.rho*vinf.Cd*(v(t))^2*vinf.Af + vinf.m*a(t));  % Torque needed at the wheels
    else
        Tw(t) = 0;
    end
    
    Ww(t) = v(t)/vinf.rwh;       % Wheel Speed Needed - [rad/sec]
    Pd(t) = Ww(t)*Tw(t);      % Power Demand       - [W]
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%---------------------Store Data into structure---------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% From the drive cycle
cyc_data.Tw = Tw;              
cyc_data.Ww = Ww;              
cyc_data.cyc_spd = cyc_mph(:,2);
cyc_data.time_cyc = time_cyc;
cyc_data.cyc_name = cyc_name;
cyc_data.Pd = Pd;
cyc_data.dt = dt;

% plot(cyc_time,a)
% xlabel('time (s)')
% ylabel('acceleration (m/s^2)')
% figure(222);
% [AX,H1,H2] = plotyy(1:time_cyc,v/mph_mps,1:time_cyc,a);
% set(get(AX(1),'Ylabel'),'String','Vehicle Speed (MPH)','fontWeight','bold','fontSize',12)
% set(get(AX(2),'Ylabel'),'String','acceleration (m/s^2)','fontWeight','bold','fontSize',12)
% set(H1,'LineWidth',2);
% set(H2,'LineWidth',2);
% set(AX(2),'fontWeight','bold','fontSize',12)
% xlabel('time (sec)','fontWeight','bold','fontSize',12);
% set(gca,'fontSize',12,'fontWeight','bold'),grid
% % % %%
% clf;
clear fig;
figure(122); clf
subplot(2,1,1)
[AX,H1,H2] = plotyy(1:time_cyc,Tw,1:time_cyc,Ww);
set(get(AX(1),'Ylabel'),'String','Wheel Torque (Nm)','fontWeight','bold','fontSize',16)
set(get(AX(2),'Ylabel'),'String','Wheel Speed (rad/sec)','fontWeight','bold','fontSize',16)
set(H1,'LineWidth',2);
set(H2,'LineWidth',2);
set(AX(2),'fontWeight','bold','fontSize',16)
set(AX,'XTick',[])
xlabel({'(a)'},'fontWeight','bold','fontSize',16);
set(gca,'fontSize',12,'fontWeight','bold'),grid
% 
% clf;
% figure(123);
% [AX,H1,H2] = plotyy(1:time_cyc,Pd,1:time_cyc,Ww);
% set(get(AX(1),'Ylabel'),'String','Power Demand','fontWeight','bold','fontSize',12)
% set(get(AX(2),'Ylabel'),'String','Wheel Speed (rad/sec)','fontWeight','bold','fontSize',12)
% set(H1,'LineWidth',2);
% set(H2,'LineWidth',2);
% set(AX(2),'fontWeight','bold','fontSize',12)
% xlabel('time (sec)','fontWeight','bold','fontSize',12);
% set(gca,'fontSize',12,'fontWeight','bold'),grid
% 
hold on
% figure(122);
subplot(2,1,2)
plot(1:time_cyc,cyc_mph(:,2),'LineWidth',2)
ylabel('Speed (mph)','fontWeight','bold','fontSize',16)
xlabel({'time (sec)','(b)'},'fontWeight','bold','fontSize',16);
set(gca,'fontSize',12,'fontWeight','bold'),grid
% title({cyc_name},'fontWeight','bold','fontSize',16)