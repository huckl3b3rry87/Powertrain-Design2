%%    
figure(17); clf;
[C,h] = contour(m_map_spd*rads2rpm, m_map_trq, m_eff_map');
hold on; 
plot(W_mot*rads2rpm,T_mot, 'ko', 'markersize', 8, 'markerf', 'g')
hold on;
plot(m_max_spd*rads2rpm, m_max_trq,  'ro', 'markersize', 12, 'markerf', 'r')
hold on
plot(m_max_spd*rads2rpm, m_max_gen_trq,  'ro', 'markersize', 12, 'markerf', 'r')
xlabel('Motor Speed (RPM)');
ylabel('Torque (Nm)');
set(gca,'FontSize',20,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',25,'fontWeight','bold')
legend('Motor Efficiency','Motor Opperating Points','Maximum Torque'),grid 
hold on
clabel(C,h, 'fontsize', 25,'fontWeight','bold');    
hold off    