function [cyc_data] = Drive_Cycle(param, vinf, cyc_name )
dt = 1;

cd('Drive_Cycle')
switch cyc_name
    % ---- Standard Cycles -=---
    case 'HWFET'  % user selects HWFET
        load CYC_HWFET;
    case 'UDDS' 
        load CYC_UDDS;
    case 'US06' 
        load CYC_US06;
    case 'LA92'  
        load CYC_LA92;
    case 'FTP75'     % Get this cycle
        load FTP75;
    case 'COMMUTER'
        load CYC_COMMUTER;
    case 'SHORT_CYC_HWFET'
        load SHORT_CYC_HWFET;
        
        % City Cycles
    case'INDIA_URBAN'
        load CYC_INDIA_URBAN_SAMPLE;
    case 'MANHATTAN';
        load CYC_MANHATTAN;
    case 'Nuremberg'
        load CYC_NurembergR36;
    case 'NYCC'
        load CYC_'NYCC';
    case 'D10'
        load D10;
        
        % ----  AV Cycles -------
    case 'US06_AV'
        load CYC_US06_AV;
    case 'HWFET_AV'
        load CYC_HWFET_AV;
    case 'D10AV'
        load D10AV;
end
cd ..

v = cyc_mph(:,2)*param.mph_mps;  % Cycle Speed in (m/s)
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
cyc_data.Tw = Tw;              
cyc_data.Ww = Ww;              
cyc_data.cyc_spd = cyc_mph(:,2);
cyc_data.time_cyc = time_cyc;
cyc_data.cyc_name = cyc_name;
cyc_data.Pd = Pd;
cyc_data.dt = dt;
cyc_data.cyc_time = cyc_time;

% --------------Can bring outside the function later to plot
% plot(cyc_time,a)
% xlabel('time (s)')
% ylabel('acceleration (m/s^2)')
% figure;
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
% figure;
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
% figure;
% [AX,H1,H2] = plotyy(1:time_cyc,Pd,1:time_cyc,Ww);
% set(get(AX(1),'Ylabel'),'String','Power Demand','fontWeight','bold','fontSize',12)
% set(get(AX(2),'Ylabel'),'String','Wheel Speed (rad/sec)','fontWeight','bold','fontSize',12)
% set(H1,'LineWidth',2);
% set(H2,'LineWidth',2);
% set(AX(2),'fontWeight','bold','fontSize',12)
% xlabel('time (sec)','fontWeight','bold','fontSize',12);
% set(gca,'fontSize',12,'fontWeight','bold'),grid
% 
% figure;
% plot(1:time_cyc,cyc_mph(:,2),'LineWidth',2)
% ylabel('Speed (mph)','fontWeight','bold','fontSize',12)
% xlabel('time (sec)','fontWeight','bold','fontSize',12);
% set(gca,'fontSize',12,'fontWeight','bold'),grid
% title('HWFET cycle','fontWeight','bold','fontSize',16)

end

