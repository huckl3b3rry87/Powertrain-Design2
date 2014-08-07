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
        Tw(t) = rwh*((m*g*atan(sin(alpha)) + Frr*m*g*atan(cos(alpha)) + 0.5*rho*Cd*(v(t))^2)*Af + (m)*a(t));  % Torque needed at the wheels
    else
        Tw(t) = 0;
    end
    
    Ww(t) = (v(t)/rwh);       % Wheel Speed Needed - [rad/sec]
    Pd(t) = Ww(t)*Tw(t);      % Power Demand       - [W]
end

% 
% % %%
% clf;
% figure(122);
% [AX,H1,H2] = plotyy(1:time_cyc,Tw,1:time_cyc,Ww);
% set(get(AX(1),'Ylabel'),'String','Wheel Torque (Nm)','fontWeight','bold','fontSize',12)
% set(get(AX(2),'Ylabel'),'String','Wheel Speed (rad/sec)','fontWeight','bold','fontSize',12)
% set(H1,'LineWidth',2);
% set(H2,'LineWidth',2);
% set(AX(2),'fontWeight','bold','fontSize',12)
% xlabel('time (sec)','fontWeight','bold','fontSize',12);
% set(gca,'fontSize',12,'fontWeight','bold'),grid
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
% figure(2);
% plot(1:time_cyc,cyc_mph(:,2),'LineWidth',2)
% ylabel('Speed (mph)','fontWeight','bold','fontSize',12)
% xlabel('time (sec)','fontWeight','bold','fontSize',12);
% set(gca,'fontSize',12,'fontWeight','bold'),grid
% title('HWFET cycle','fontWeight','bold','fontSize',16)