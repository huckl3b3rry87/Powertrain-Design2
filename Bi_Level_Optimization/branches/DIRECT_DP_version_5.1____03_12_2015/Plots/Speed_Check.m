%%  Motor is always connected to vehicle...
figure(14);clf;

We = W_eng./(GEAR_save*FD);
Wm = W_mot./(FD*G);

plot(cyc_time,Ww,'color','r','LineWidth',8);
hold on 
plot(cyc_time,Wm,'k--','LineWidth',4);
hold on
plot(cyc_time,We,'color','g','LineWidth',2);
hold on
legend('Wheel','Motor (at road)','Engine (at road)');
ylabel('Speed (rad/s)','fontWeight','bold','fontSize',20)
xlabel('time (sec)','fontWeight','bold','fontSize',20);
set(gca,'fontSize',20,'fontWeight','bold'),grid
hold off