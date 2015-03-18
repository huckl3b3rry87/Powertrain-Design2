%%
n = 4;
q = 1;
h = 9;
y = 12;
t = 1.25*cyc_data.time_cyc;
figure(34); clf;
subplot(4,1,1);
GRAPH =  plot(cyc_data.cyc_time,sim.CO,'rx-','LineWidth',2,'markersize', 3); grid on;
set(GRAPH,'marker','x', 'markersize', n, 'markerf', 'b','linewidth',q);
% H = legend('Engine Speed','Motor Speed');set(H,'fontsize', h)
ylabel({'CO rate', '(g/s)'},'fontWeight','bold','fontSize',y)
set(gca,'fontSize',y,'fontWeight','bold');   
% set(gca, 'pos', [0.129559619306594 0.425027519818801 0.775 0.172933975084938]);
set(gca,'XTickLabel',[]);
xlim([0 t]);

subplot(4,1,2);
GRAPH =  plot(cyc_data.cyc_time,sim.NOx,'rx-','LineWidth',2,'markersize', 3); grid on;
set(GRAPH,'marker','x', 'markersize', n, 'markerf', 'b','linewidth',q);
% H = legend('Engine Speed','Motor Speed');set(H,'fontsize', h)
ylabel({'NOx rate', '(g/s)'},'fontWeight','bold','fontSize',y)
set(gca,'fontSize',y,'fontWeight','bold');   
% set(gca, 'pos', [0.129559619306594 0.425027519818801 0.775 0.172933975084938]);
set(gca,'XTickLabel',[]);
xlim([0 t]);

subplot(4,1,3);
GRAPH =  plot(cyc_data.cyc_time,sim.HC,'rx-','LineWidth',2,'markersize', 3); grid on;
set(GRAPH,'marker','x', 'markersize', n, 'markerf', 'b','linewidth',q);
% H = legend('Engine Speed','Motor Speed');set(H,'fontsize', h)
ylabel({'HC rate', '(g/s)'},'fontWeight','bold','fontSize',y)
set(gca,'fontSize',y,'fontWeight','bold');   
% set(gca, 'pos', [0.129559619306594 0.425027519818801 0.775 0.172933975084938]);
set(gca,'XTickLabel',[]);
xlim([0 t]);
% TOP = [max(sim.W_eng*rads2rpm), max(sim.W_mot*rads2rpm)];
% ylim([0 1.1*max(TOP)])

subplot(4,1,4);
GRAPH =  plot(cyc_data.cyc_time,sim.inst_fuel,'rx-','LineWidth',2,'markersize', 3); grid on;
set(GRAPH,'marker','x', 'markersize', n, 'markerf', 'b','linewidth',q);
% H = legend('Engine Speed','Motor Speed');set(H,'fontsize', h)
ylabel({'Fuel rate', '(g/s)'},'fontWeight','bold','fontSize',y)
set(gca,'fontSize',y,'fontWeight','bold');   
% set(gca, 'pos', [0.129559619306594 0.425027519818801 0.775 0.172933975084938]);
set(gca,'XTickLabel',[]);
xlim([0 t]);

xlabel('time (sec)');
