% Natralisitic Drive Cycle
%%
clear
%%

for i = 1:1:length(VarName1)
    cyc_mph(i,1) = i;  % Time
    cyc_mph(i,2) = VarName1(i); % Speed
end

time_cyc =  length(cyc_mph(:,1));
figure(2);
plot(1:time_cyc,cyc_mph(:,2),'LineWidth',2)
ylabel('Speed (mph)','fontWeight','bold','fontSize',12)
xlabel('time (sec)','fontWeight','bold','fontSize',12);
set(gca,'fontSize',12,'fontWeight','bold'),grid
title('Ann Arbor cycle','fontWeight','bold','fontSize',16)

save('D10','cyc_mph')