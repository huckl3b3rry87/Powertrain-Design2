%%  Motor is always connected to vehicle...
figure(14);clf;

We = sim.W_eng./(vinf.gear(sim.GEAR)*dvar.FD);
Wm = sim.W_mot./(dvar.FD*dvar.G);
Ww = cyc_data.Ww;
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