Te = T_eng.*GEAR_save*FD;
Tm = T_mot.*GEAR_save*FD*G;
T_total = Te + Tm;

figure(13);clf;
plot(cyc_time,Te,'color','g','LineWidth',2);
hold on 
plot(cyc_time,Tm,'LineWidth',3);
hold on
plot(cyc_time,Tw,'color','y','LineWidth',6);
hold on
plot(cyc_time,T_total,'r--','LineWidth',2)
legend('Engine (at road)','Motor (at road)','T required','Motor + Engine')
ylabel('Torque (Nm)','fontWeight','bold','fontSize',20);grid
xlabel('time (sec)','fontWeight','bold','fontSize',20);
set(gca,'fontSize',20,'fontWeight','bold')
hold off    