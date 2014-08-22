figure(10);clf
subplot(3,1,1), plot(time_sim,x_sim)
xlabel('Time (s)')
ylabel('Distance (m)')
grid
   
subplot(3,1,2), plot(time_sim,V_sim)
xlabel('Time (s)')
ylabel('Velocity (MPH)')
grid
   
subplot(3,1,3), plot(time_sim,acc_sim)
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
grid
   