
h = 10;
r = 12;
l = 1.5;
figure(10);clf
subplot(3,1,1), plot(time_sim,x_sim,'k-','linewidth',l)
ylabel({'Distance', '(m)'})
set(gca,'FontSize',h,'fontWeight','bold')
set(gca,'XTickLabel',[]);
set(findall(gcf,'type','text'),'FontSize',r,'fontWeight','bold')
grid
   
subplot(3,1,2), plot(time_sim,V_sim,'k-','linewidth',l)
ylabel({'Velocity','(mph)'})
set(gca,'FontSize',h,'fontWeight','bold')
set(gca,'XTickLabel',[]);
set(findall(gcf,'type','text'),'FontSize',r,'fontWeight','bold')
grid
   
subplot(3,1,3), plot(time_sim,acc_sim,'k-','linewidth',l)
xlabel('Time (s)')
ylabel({'Acceleration', '(m/s^2)'})
set(gca,'FontSize',h,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',r,'fontWeight','bold')
grid
   