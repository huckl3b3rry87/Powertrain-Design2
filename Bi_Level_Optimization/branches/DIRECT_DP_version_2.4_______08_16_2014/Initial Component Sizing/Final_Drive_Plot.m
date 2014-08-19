figure(2);clf
plot(V_6/param.mph_mps,FD_sim,'k.','markersize',25)
hold on
plot(V_act/param.mph_mps,FD_sim,'g*')
hold on
plot(Max_Spd, FD_sim(I), 'ro','markersize',30)
hold on
legend('Gear Six', 'All Gears', 'Selected FD Ratio')
xlabel('Vehicle Speed (MPH)');
ylabel('Final Drive Ratio');
title(['Final Drive Ratio Selection = ', num2str(FD_sim(I))])
set(gca,'FontSize',20,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',25,'fontWeight','bold')
grid; hold off