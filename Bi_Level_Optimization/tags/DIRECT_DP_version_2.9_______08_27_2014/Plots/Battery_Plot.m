%%    
figure(30); clf;
plot(vinf.ess_soc,vinf.ess_max_pwr_dis/1000, 'b','linewidth', 12)
hold on;
plot(vinf.ess_soc, -vinf.ess_max_pwr_chg/1000, 'r', 'linewidth', 12)
hold on
plot( sim.SOC_final, sim.Pbatt_sim/1000,   'ko', 'markersize',10, 'markerf', 'g','linewidth',3)
xlabel('SOC');
ylabel('Power (kW)');
set(gca,'FontSize',20,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',25,'fontWeight','bold')
legend('Maximum Battery Discharge Power','Maximum Battery Charge Power','Opperating Points'),grid 


% figure(34); clf
% plot(cyc_time,i)
% legend('current')
% 
% figure(35); clf
% plot(cyc_time,C_rate)
% legend('C Rate')