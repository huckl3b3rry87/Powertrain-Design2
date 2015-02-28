
for t = 1:length(cyc_time)
    if Pd(t) < 0 
        R(t) = 100;
    else
        R(t) = 0;
    end
end

figure(33); clf
plot(cyc_time,T_mot,'r','linewidth',3)
hold on
plot(cyc_time,Tm_actual_sim,'k')
hold on
plot(cyc_time, T_brake_sim,'g','linewidth',3),grid
hold on
plot(cyc_time,R)
hold on 
plot(cyc_time,ENG_state*150,'m','linewidth',3)
legend('Motor Torque Reqested','Motor Torque Actual','Brake Torque','Pd < 0','ENG state')