figure(16); clf;
[C,h] = contourf(eng_consum_spd*rads2rpm, eng_consum_trq, eng_bsfc',  [216 300:10:350, 275, 350:25:800]);
clabel(C,h, 'fontsize', 8);
axis([0 5500 0 max(eng_max_trq+5)]);
hold on;
plot(eng_map_spd*rads2rpm, eng_max_trq, 'r', 'linewidth', 4);
hold on 
plot(W_eng*rads2rpm,T_eng, 'ko', 'markersize', 12, 'markerf', 'g'),grid
xlabel('Engine Speed (RPM)');
ylabel('Torque (Nm)');
set(gca,'FontSize',20,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',25,'fontWeight','bold')
legend('BSFC (g/kWh)','Maximum Engine Torque','Engine Opperating Points')    
    
