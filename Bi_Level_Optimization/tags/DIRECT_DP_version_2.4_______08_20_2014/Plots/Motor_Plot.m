%%    
eff = linspace(min(min(mc_eff_map)),max(max(mc_eff_map)),15);

figure(17); clf;
[C,h] = contourf(m_map_spd*rads2rpm, m_map_trq, m_eff_map',eff);
hold on; 
plot(W_mot*rads2rpm,Tm_actual_sim, 'ko', 'markersize', 8, 'markerf', 'g','linewidth',3)
hold on;
plot(m_map_spd *rads2rpm, m_max_trq,  'k','linewidth',12)
hold on
plot(m_map_spd *rads2rpm, m_max_gen_trq, 'k','linewidth',12)
xlabel('Motor Speed (RPM)');
ylabel('Torque (Nm)');
set(gca,'FontSize',20,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',25,'fontWeight','bold')
legend('Motor Efficiency','Motor Opperating Points','Maximum Continuous Torque'),grid 

colorbar('FontSize',17,'fontWeight','bold')

% colorbar('YTickLabel',{num2str(eff(1)),num2str(eff(2)),num2str(eff(3)),num2str(eff(4)),num2str(eff(5)),num2str(eff(6))...
%     num2str(eff(7)),num2str(eff(8)),num2str(eff(9)),num2str(eff(10)), num2str(eff(11)),num2str(eff(12)),num2str(eff(13)),num2str(eff(14)),num2str(eff(15))},'FontSize',17,'fontWeight','bold','YTick',eff)

hold off    