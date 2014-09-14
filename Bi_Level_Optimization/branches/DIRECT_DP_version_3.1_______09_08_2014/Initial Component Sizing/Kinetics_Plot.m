
h = 14;
r = 16;
figure(11);clf
subplot(2,2,1), plot(time_sim,Teng_sim,'k-','linewidth',3) 
xlabel('Time (s)')
ylabel({'Engine Torque', '(Nm)'})
set(gca,'FontSize',h,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',r,'fontWeight','bold')
grid

subplot(2,2,2), plot(time_sim,Tm_sim,'k-','linewidth',3) 
xlabel('Time (s)')
ylabel({'Motor Torque', '(Nm)'})
set(gca,'FontSize',h,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',r,'fontWeight','bold')
grid

subplot(2,2,3), plot(time_sim,We_sim,'k-','linewidth',3) 
xlabel('Time (s)')
ylabel({'Engine Speed', '(rpm)'})
set(gca,'FontSize',h,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',r,'fontWeight','bold')
grid

subplot(2,2,4), plot(time_sim,Wm_sim,'k-','linewidth',3) 
xlabel('Time (s)')
ylabel({'Motor Speed', '(rpm)'})
set(gca,'FontSize',h,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',r,'fontWeight','bold')
grid
