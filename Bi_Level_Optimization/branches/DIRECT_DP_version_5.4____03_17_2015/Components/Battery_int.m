%% Battery Parameters
%% Chiao-Ting 10/10/2013
%% 
%% Reference: 
%% T. Gray and M. Shirk, "2010 Toyota Prius VIN 0462 Hybrid Electric 
%% Vehicle Battery Test Results," The Idaho National Laboratory, Report 
%% #INL/EXT-13-28025, 2013.
%%
%% battery of 3rd generation Prius (2010 Prius): 6.5 Ah, 1.3~1.35 kWh
%% ================================================================== %%

% Battery modules do not affect this battery as is!!!

ess_cap_ah = 6.5; % (A*h)

ess_wh = 0:25:1350; % [wh]
ess_soc = ess_wh/1350; % [-], dsoc ~ 0.0185

%% whole battery pack!! not single module
ess_r_dis = [0.430 	0.428 	0.427 	0.425 	0.423 	0.422 	0.420 	0.418 	0.417 	0.415 	0.413 	0.412 	0.410 	0.411 	0.412 	0.413 	0.413 	0.414 	0.415 	0.416 	0.417 	0.418 	0.418 	0.419 	0.420 	0.421 	0.422 	0.423 	0.424 	0.425 	0.426 	0.427 	0.428 	0.431 	0.435 	0.438 	0.442 	0.445 	0.448 	0.452 	0.455 	0.461 	0.466 	0.472 	0.478 	0.483 	0.489 	0.494 	0.500 	0.515 	0.529 	0.544 	0.558 	0.573 	0.587 ];
ess_r_chg = [0.470 	0.465 	0.460 	0.455 	0.450 	0.445 	0.440 	0.438 	0.437 	0.435 	0.433 	0.432 	0.430 	0.430 	0.430 	0.430 	0.430 	0.430 	0.431 	0.433 	0.434 	0.436 	0.437 	0.439 	0.440 	0.442 	0.443 	0.445 	0.446 	0.448 	0.449 	0.451 	0.452 	0.455 	0.458 	0.461 	0.464 	0.467 	0.473 	0.480 	0.486 	0.492 	0.499 	0.505 	0.513 	0.521 	0.530 	0.538 	0.546 	0.559 	0.578 	0.596 	0.614 	0.632 	0.650 ];
ess_voc = [170.000 	187.000 	195.000 	199.250 	201.640 	203.400 	204.280 	205.510 	206.470 	207.460 	208.410 	209.190 	209.910 	210.660 	211.290 	211.820 	212.390 	212.900 	213.350 	213.770 	214.150 	214.470 	214.770 	215.030 	215.290 	215.530 	215.750 	215.940 	216.120 	216.270 	216.420 	216.570 	216.710 	216.860 	217.000 	217.160 	217.290 	217.390 	217.560 	217.760 	217.960 	218.260 	218.660 	219.020 	219.360 	219.860 	220.460 	221.080 	221.820 	222.820 	224.020 	225.420 	227.100 	229.000 	233.020];

% % For Battery sizing
% Rint_size = [ess_r_dis(4)/module_number];
% Voc_size =  [ess_voc(4)/module_number];

%% the resistance down scaling helps to match the 100kW net power spec of 2010 Prius (mentioned in Namwook Kim's journal)
ess_max_pwr_dis = (ess_voc.^2)./(4*ess_r_dis)*0.98; % 110~105kW in 40-70 percent SOC window
ess_max_pwr_chg = (ess_voc.^2)./(4*ess_r_chg)*0.98; % 105~100kW in 40-70 percent SOC window

% plot(ess_soc,ess_max_pwr_dis/1000)