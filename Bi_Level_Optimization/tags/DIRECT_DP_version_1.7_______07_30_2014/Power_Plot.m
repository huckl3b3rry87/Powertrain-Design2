
We = W_eng./(GEAR_save*FD);
Wm = W_mot./(GEAR_save*FD*G);

% plot(cyc_time,Ww,'color','r','LineWidth',8);
% hold on 
% plot(cyc_time,Wm,'k--','LineWidth',4);
% hold on
% plot(cyc_time,We,'color','g','LineWidth',2);
% hold on
% legend('Wheel','Motor (at road)','Engine (at road)');
% ylabel('Speed (rad/s)','fontWeight','bold','fontSize',20)
% xlabel('time (sec)','fontWeight','bold','fontSize',20);
% set(gca,'fontSize',20,'fontWeight','bold'),grid
% hold off


Te = T_eng.*GEAR_save*FD;
Tm = T_mot.*GEAR_save*FD*G;
T_total = Te + Tm;


Peng =  Te.*We;
Pbatt_sim = Pbat;
Pd_tot = Pd + Paux;




% figure(13);clf;
% plot(cyc_time,Te,'color','g','LineWidth',2);
% hold on 
% plot(cyc_time,Tm,'LineWidth',3);
% hold on
% plot(cyc_time,Tw,'color','y','LineWidth',6);
% hold on
% plot(cyc_time,T_total,'r--','LineWidth',2)
% legend('Engine (at road)','Motor (at road)','T required','Motor + Engine')
% ylabel('Torque (Nm)','fontWeight','bold','fontSize',20);grid
% xlabel('time (sec)','fontWeight','bold','fontSize',20);
% set(gca,'fontSize',20,'fontWeight','bold')
% hold off    