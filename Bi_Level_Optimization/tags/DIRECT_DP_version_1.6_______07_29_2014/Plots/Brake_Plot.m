figure(33); clf
plot(cyc_time,T_mot,'r','linewidth',3)
hold on
plot(cyc_time,Tm_actual_sim,'k')
hold on
plot(cyc_time, T_brake_sim,'g','linewidth',3),grid
legend('Motor Torque Reqested','Motor Torque Actual','Brake Torque')