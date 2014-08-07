%% Engine Parameters
%% Xiaowu 2011/05/13
%% 
%% 2010 Prius (third generation Prius)
%% from official website:
%%  power: 98 hp @ 5200 rpm (73 kW @ 5200 rpm)
%%  torque: 105 lb.-ft. @ 4000 rpm (142 N·m @ 4000 rpm)
%%  hybrid system net power: 134 hp (100 kW)
%% ================================================================== %%

rpm2rads = pi/30;
rads2rpm = 1/rpm2rads;
mph2mps = 1609/3600;
mps2mph = 1/mph2mps;
g2gallon=1/841.4/3.785;


e_inertia=0.18; % [kg*m^2], rotational inertia of the engine (unknown)
eng_consum_spd = [0,1000:200:5200]*rpm2rads;
eng_consum_trq=[-6,0:10:150];
eng_bsfc= [
%-6    0    10   20    30    40   50   60   70  80   90   100  110  120  130   140  150, [Nm]
 0    800  700  590  425   278  260  250  240  235  230  228  235  245  249   255  260;...%0, [rpm]
 0    800  700  592  424   280  254  243  230  225  221  226  230  236  243   245  247;...%1000
 0    800  700  594  423   281  254  243  230  224  219  223  228  232  241   243  246;...%1200
 0    800  700  596  422   282  253  242  229  224  217  221  225  230  238   242  245;...%1400
 0    800  680  600  421   283  253  242  229  223  216  219  223  228  235   240  244;...%1600
 0    800  680  600  420   284  253  242  230  223  218  217  221  226  230   238  242;...%1800
 0    800  690  600  420   285  254  242  231  224  219  215  219  226  228   236  241;...%2000
 0    800  690  600  420   286  254  243  232  225  220  216  218  225  227   234  238;...%2200
 0    800  690  600  420   287  255  244  233  226  221  217  217  224  226   232  234;...%2400
 0    800  700  600  421   288  256  244  234  227  222  218  218  224  226   230  232;...%2600
 0    800  700  610  424   289  257  245  236  228  223  219  219  225  227   230  233;...%2800
 0    800  700  612  428   289  258  246  238  229  224  221  220  226  228   232  235;...%3000
 0    800  700  615  430   289  259  246  240  231  226  222  222  227  228   234  238;...%3200
 0    800  700  617  433   292  261  247  242  233  227  224  224  228  229   235  240;...%3400
 0    800  700  620  437   295  263  248  245  235  229  226  227  228  231   236  242;...%3600
 0    800  700  622  439   300  265  249  247  237  230  227  229  229  233   237  244;...%3800
 0    800  700  624  441   305  267  250  249  239  233  230  231  232  235   239  246;...%4000
 0    800  700  628  446   310  270  252  250  241  237  231  233  235  237   241  248;...%4200
 0    800  700  631  448   315  272  254  251  244  241  235  236  238  239   242  249;...%4400
 0    800  700  635  450   320  274  255  252  248  242  242  239  240  241   243  251;...%4600
 0    800  710  638  452   325  276  257  254  250  244  244  241  242  242   244  253;...%4800
 0    800  710  640  454   330  279  259  255  252  246  245  243  244  243   244  255;...%5000 
 0    800  700  640  455   335  282  262  257  253  248  246  244  245  245   246  257];  %5200
 
engine_p = eng_consum_spd'*eng_consum_trq;
eng_consum_fuel = eng_bsfc.*(engine_p/1000)/3600;

eng_consum_fuel_raw = eng_consum_fuel;

eng_consum_fuel(1,:) = zeros(1,length(eng_consum_trq));
eng_consum_fuel(:,1) = zeros(length(eng_consum_spd),1);
eng_consum_fuel(:,2) = zeros(length(eng_consum_spd),1);
eng_consum_fuel(1,:) = eng_consum_fuel(2,:);
eng_consum_fuel(:,2) = 0.5*eng_consum_fuel(:,3);
eng_consum_fuel(1,2) = 0;

eng_map_spd = [0,1000:200:5200]*rpm2rads;
eng_max_trq = [85   105  108  112  116  120  123  127  131  133  135  137  138  140  142 142.5 143 141.5 140  139 138  136  135];

%% optimal BSFC
We_optbsfc = (0:50:5200)*rpm2rads;
Te_optbsfc = [0 	0 	0 	0 	0 	0 	0 	0 	0 	0 	0 	0 	0 	0 	0 	0 	0 	0 	0 	0 	86 	89 	90 	90 	90 	90 	90 	90 	90 	90 	90 	90 	90 	91 	92 	95 	97 	99 	100 	100 	100 	100 	100 	100 	100 	100 	100 	103 	107 	110 	110 	110 	110 	110 	110 	110 	110 	110 	110 	110 	110 	110 	110 	110 	110 	110 	110 	110 	111 	113 	115 	116 	119 	120 	120 	121 	123 	126 	128 	129 	130 	130 	130 	130 	130 	130 	130 	130 	130 	134 	135 	136 	137 	138 	139 	138 	138 	138 	137 	137 	136 	136 	136 	135 	135];
fuel_optbsfc = interp2(eng_consum_trq,eng_consum_spd,eng_consum_fuel,Te_optbsfc,We_optbsfc);


% %% ==============================
% figure(100); clf;
% [C,h] = contour(eng_consum_spd*rads2rpm, eng_consum_trq, eng_bsfc', [216 200:5:250, 275, 300:100:800]);
% clabel(C,h, 'fontsize', 8);
% axis([0 5500 0 150]);
% 
% hold on;
% plot(eng_map_spd*rads2rpm, eng_max_trq, 'k', 'linewidth', 2);
% plot(We_optbsfc*rads2rpm, Te_optbsfc, 'ro-', 'markersize', 3, 'markerf', 'r');
% hold on 
% plot(W_eng*rads2rpm,T_eng, 'ko', 'markersize', 12, 'markerf', 'g')
% xlabel('Engine Speed (RPM)');
% ylabel('Torque (Nm)');
% title('Engine Opperating Points');
% set(gca,'FontSize',20,'fontWeight','bold')
% set(findall(gcf,'type','text'),'FontSize',25,'fontWeight','bold')


%% ==============================
% % figure; clf;
% % mesh(eng_consum_spd, eng_consum_trq, eng_bsfc');
% % xlabel('Spd');
% % ylabel('Trq')
% % zlabel('BSFC (g/kWh)');
% % 
% % % ==============================
% % figure; clf;
% % mesh(eng_consum_spd, eng_consum_trq, eng_consum_fuel_raw');
% % xlabel('Spd');
% % ylabel('Trq')
% % zlabel('Fuel Rate (g/s)');



