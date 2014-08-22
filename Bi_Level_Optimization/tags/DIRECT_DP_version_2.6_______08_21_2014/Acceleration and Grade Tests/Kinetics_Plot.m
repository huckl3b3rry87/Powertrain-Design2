figure(11);clf
subplot(2,2,1), plot(time_sim,Teng_sim) 
xlabel('Time (s)')
ylabel('Engine Torque (Nm)')
grid

subplot(2,2,2), plot(time_sim,Tm_sim) 
xlabel('Time (s)')
ylabel('Motor Torque (Nm)')
grid

subplot(2,2,3), plot(time_sim,We_sim) 
xlabel('Time (s)')
ylabel('Engine Speed (RPM)')
grid

subplot(2,2,4), plot(time_sim,Wm_sim) 
xlabel('Time (s)')
ylabel('Motor Speed (RPM)')
grid
