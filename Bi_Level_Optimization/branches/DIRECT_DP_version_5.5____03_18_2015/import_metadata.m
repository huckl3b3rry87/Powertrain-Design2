clear

filename = 'Query1.xlsx';
sheet = 1;

driver = 'A:A';
trip = 'B:B';
time = 'C:C';
speed = 'D:D';

t_temp = xlsread(filename, sheet, time);
s_temp = xlsread(filename, sheet, speed);
driver_t = xlsread(filename, sheet, driver);
trip_t = xlsread(filename, sheet, trip);

%% Break it up
clear z set tr_n tr_o
%Need size and location of the trips
z = 1;
tr_o = 0;
set(z) = 1;
set_w(z) = 0;
for i = 1:length(t_temp)
    tr_n = trip_t(i);
    
    if  i ~=1 && tr_o ~= tr_n         % There is no change in the last trip type!
          set_w(z+1) = i - sum(set_w(1:z));
          set(z+1) = i;
        z = z+1;
    end
    tr_o = tr_n;
end

% Find the max
[B,I] = sort(set_w);

bottom = length(B)-10;

a = 0;
index(1) = [1];
speed_aa = [];
for w = 1:10
    I_start(w)=set(I(end - a)-1);
    I_end(w) = set(I(end - a)) - 1;
    a = a + 1;
end

% find the driver and the test run
for i = 1:10 
    time_sim = 1:length(s_temp(I_start(i):I_end(i)));
    speed_sim = s_temp(I_start(i):I_end(i));
    driver_save(i) = driver_t(I_start(i));
    trip_save(i) = trip_t(I_start(i));
    figure(i);
    plot(time_sim,speed_sim*2.23694,'LineWidth',2);
    ylabel('Speed (mph)','fontWeight','bold','fontSize',12)
    xlabel('time (sec)','fontWeight','bold','fontSize',12);
    title('Ann Arbor cycle','fontWeight','bold','fontSize',16)
    
    %     % Save the cycle for acceleration requirements
    %     spd_save{i} = speed_sim;
    %     time_save{i} = time_sim;

    cyc_mph(:,1) =  time_sim;
    cyc_mph(:,2) = speed_sim*2.23694;
    
    eval(sprintf( 'AA_D%d_T%d = cyc_mph', driver_save(i), trip_save(i)))
    save(sprintf( 'AA_D%d_T%d', driver_save(i), trip_save(i)),'cyc_mph')
    clear  cyc_mph
end

% figure(22)
% plot(cell2mat(time_save(1)),cell2mat(spd_save(1)))
% close all

% Modified Drive Cycle
% Refereing to EE, this cycle consists of four different parts 
%1. from 19 seconds to 390 seconds of D32 T187
%2. from 36 sec to  180 sec of D46 T62
%3. from 64 sec to 204 sec of D104 T52
V_new = [s_temp(I_start(1)+19:(I_start(1)+390));s_temp((I_start(9)+36):(I_start(9)+180));s_temp(I_start(5)+64:I_start(5)+204)];
time_new = 1:length(V_new);
figure(33);clf
plot(time_new,V_new*2.23694,'LineWidth',2)
ylabel('Speed (mph)','fontWeight','bold','fontSize',12)
xlabel('time (sec)','fontWeight','bold','fontSize',12);
title('Natralistic Ann Arbor cycle','fontWeight','bold','fontSize',16)
grid on

cyc_mph(:,1) =  time_new;
cyc_mph(:,2) = V_new*2.23694;
save('AA_final','cyc_mph');

% EE =
% 
%     32   187  1 
%     19     3  2
%      4   258  3
%     52   100  4
%    104    52  5
%     52    62  6
%      4   286  7 
%    104    56  8
%     46    62  9
%    104   195  10
