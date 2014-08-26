%%
n = 4;
q = 1;
h = 9;
y = 12;
t = 1.15*cyc_data.time_cyc;
figure(1); clf;
subplot(5,1,1);
[AX,L1,L2] = plotyy(cyc_data.cyc_time, cyc_data.cyc_spd, cyc_data.cyc_time,sim.SOC_final);
set(L1,'marker','x', 'markersize', n, 'markerf','b','linewidth',q);
set(L2,'marker','x', 'markersize', n, 'markerf','g','linewidth',q);
set(AX, 'fontsize', y,'fontweight','bold','XTickLabel',[]);
ylabel({'Speed', '(mph)'}, 'parent', AX(1));
ylabel('SOC', 'parent', AX(2));grid
set(AX(1),'XLim',[0 t]);
set(AX(2),'XLim',[0 t]);
set(AX(1),'YLim',[0 (max(cyc_data.cyc_spd + 5))]);
set(AX(2),'YLim',[0.39 0.82]);
H = legend('Drive Profile','SOC'); set(H,'fontsize', h)
set(gcf, 'units', 'inch', 'pos', [0.13 0.778029445073613 0.775 0.146970554926388]);

subplot(5,1,2);
[AX,L1,L2] = plotyy(cyc_data.cyc_time, sim.inst_fuel,cyc_data.cyc_time,sim.ENG_state);
set(L1,'marker','x', 'markersize', n, 'markerf','b','linewidth',q);
set(L2,'marker','x', 'markersize', n, 'markerf','g','linewidth',q);
set(AX, 'fontsize', y,'fontweight','bold','XTickLabel',[]);
ylabel({'Fuel Rate', '( g/s)'}, 'parent', AX(1));
ylabel({'Engine State','(on/off)'}, 'parent', AX(2));grid
set(AX(1),'XLim',[0 t]);
set(AX(2),'XLim',[0 t]);
set(AX(1),'YLim',[-0.1 1.2*max(sim.inst_fuel)]);
set(AX(2),'YLim',[-0.1 1.1]);
H = legend('Fuel Rate','Engine State');set(H,'fontsize', h)
set(gcf, 'units', 'inch', 'pos', [0.13 0.61328593776993 0.775 0.124322033898305]);

subplot(5,1,3);
GRAPH =  plot(cyc_data.cyc_time,sim.W_eng*rads2rpm, 'kx-',cyc_data.cyc_time,sim.W_mot*rads2rpm, 'rx-','LineWidth',2,'markersize', 3); grid on;
set(GRAPH,'marker','x', 'markersize', n, 'markerf', 'b','linewidth',q);
H = legend('Engine Speed','Motor Speed');set(H,'fontsize', h)
ylabel({'Speed', '(RPM)'},'fontWeight','bold','fontSize',y)
set(gca,'fontSize',y,'fontWeight','bold');   
set(gca, 'pos', [0.129559619306594 0.425027519818801 0.775 0.172933975084938]);
set(gca,'XTickLabel',[]);
xlim([0 t]);
TOP = [max(sim.W_eng*rads2rpm), max(sim.W_mot*rads2rpm)];
ylim([0 1.1*max(TOP)])

subplot(5,1,4);
GRAPH =  plot(cyc_data.cyc_time,sim.T_brake_sim,'m',cyc_data.cyc_time,sim.Tm_actual_sim,'g',cyc_data.cyc_time,sim.T_eng,'b'); grid on;
set(GRAPH,'marker','x', 'markersize', n, 'markerf', 'b','linewidth',q);
H=legend('Brake Torque','Motor Torque','Engine Torque');set(H,'fontsize', h)
ylabel({'Torque' ,'(Nm)'},'fontWeight','bold','fontSize',y)
set(gca,'fontSize',y,'fontWeight','bold');   
set(gca, 'pos',[0.129559619306594 0.255522536806342 0.775 0.1495]);
set(gca,'XTickLabel',[]);
END = [max(sim.T_mot), max(sim.T_eng)];
ylim([(min(sim.T_mot)-5) 1.1*(max(END))]);
xlim([0 t]);

subplot(5,1,5);
[AX,L1,L2] = plotyy(cyc_data.cyc_time, cyc_data.cyc_spd, cyc_data.cyc_time, sim.GEAR);
set(L1,'marker','x', 'markersize', n, 'markerf', 'b','linewidth',q);
set(L2,'marker','x', 'markersize', n, 'markerf', 'g','linewidth',q);
set(AX, 'fontsize', y,'fontweight','bold');
ylabel({'Speed','(mph)'}, 'parent', AX(1));
ylabel('Gear', 'parent', AX(2));grid; 
set(AX(1),'XLim',[0 t]);
set(AX(2),'XLim',[0 t]);
set(AX(1),'YLim',[0 (max(cyc_data.cyc_spd) + 5)]);
set(AX(2),'YLim',[0.5 6.5]);
xlabel('time (sec)');
H = legend('Drive Profile','Optimal Gear'); set(H,'fontsize', h)   