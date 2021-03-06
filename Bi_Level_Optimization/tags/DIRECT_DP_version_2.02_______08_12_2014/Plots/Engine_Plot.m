bsfc_temp = eng_bsfc(eng_bsfc~=0);
num1_2 = 8;
num2_3 = 6;
BSFC_L1 = linspace(min(bsfc_temp),1.09*min(bsfc_temp),8);
BSFC_L2_temp = linspace(1.09*min(bsfc_temp),1.5*min(bsfc_temp),num1_2);
BSFC_L2  = BSFC_L2_temp(2:num1_2);
BSFC_U_temp  = linspace(1.5*min(bsfc_temp),max(max(bsfc_temp),1.9*min(bsfc_temp)),num2_3);
BSFC_U  = BSFC_U_temp(2:num2_3);

BSFC = [BSFC_L1, BSFC_L2, BSFC_U];

figure(16); clf;
[C,h] = contourf(eng_consum_spd*rads2rpm, eng_consum_trq, eng_bsfc',BSFC);
axis([0 5500 0 max(eng_max_trq+5)]);
hold on;
plot(eng_map_spd*rads2rpm, eng_max_trq, 'k', 'linewidth', 8);
hold on 
plot(W_eng*rads2rpm,T_eng, 'mo', 'markersize', 10, 'markerf', 'k','linewidth',2),grid
xlabel('Speed (RPM)');
ylabel('Torque (Nm)');
title('Engine')
set(gca,'FontSize',16,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',22,'fontWeight','bold')
legend('BSFC (g/kWh)','Maximum Torque','Opperating Points')   

h = gcf;
load('EngineColorMap','mycmap')
set(h,'Colormap',mycmap)
colorbar('FontSize',16,'fontWeight','bold','YDir','reverse')

% set(h,'Colormap',flipud(mycmap))
% colorbar('FontSize',16,'fontWeight','bold',)
%% To save
% h = gcf;
% mycmap = get(h,'Colormap');
% save('EngineColorMap','mycmap');
