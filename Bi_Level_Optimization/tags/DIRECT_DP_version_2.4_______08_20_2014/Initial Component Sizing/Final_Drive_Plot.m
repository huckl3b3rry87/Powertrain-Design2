figure(2);clf

plot(V_6/param.mph_mps,FD_sim,'k.','markersize',25)
hold on
plot(V_act/param.mph_mps,FD_sim,'g*')
hold on
plot(Max_Spd_temp/param.mph_mps, FD_sim(I_temp), 'b.','markersize',40)
hold on
plot(Max_Spd/param.mph_mps, FD_sim(I), 'r.','markersize',40)
hold on
legend('Gear Six', 'All Gears','n_h', 'Selected FD Ratio')
xlabel('Vehicle Speed (MPH)');
ylabel('Final Drive Ratio');
title(['Final Drive Ratio Selection = ', num2str(FD_sim(I_temp)),'*0.9 ','= ',num2str(FD_final)])
set(gca,'FontSize',20,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',25,'fontWeight','bold')
grid; hold off