% cyc_name = 'HWFET';
% cyc_name = 'UDDS';
% cyc_name = 'US06';
% cyc_name = 'SHORT_CYC_HWFET';
% cyc_name = 'RAMP';
% cyc_name = 'LA92';
% cyc_name = 'CONST_65';
% cyc_name = 'CONST_45';
% cyc_name = 'COMMUTER';

% City
for q = 1:4
    if q == 1
        cyc_name = 'INDIA_URBAN';
    elseif q == 2
        cyc_name = 'MANHATTAN';
    elseif q == 3
        cyc_name = 'Nuremberg';
    else
        cyc_name = 'NYCC';
    end
    
    [cyc_data] = Drive_Cycle(param, vinf, cyc_name )
    
    Pd_no_regen(cyc_data.Pd > 0) = cyc_data.Pd(cyc_data.Pd > 0);
    Pd_ave = sum(Pd_no_regen)/length(Pd_no_regen)*ones(size(cyc_data.Pd));
    
    figure(q+12);clf
    plot(1:cyc_data.time_cyc,cyc_data.Pd/1000,'k','LineWidth',3)
    hold on
    plot(1:cyc_data.time_cyc,Pd_ave/1000,'b','LineWidth',3);
    legend('Instantaneous Power','Average Power')
    title(['Average Power with Zero',10,'Regenerative Braking = ',num2str(round(Pd_ave(1)/1000)),' (kW)',10,cyc_name],'fontWeight','bold','fontSize',15);
    ylabel('Power (kW)')
    xlabel('time (sec)','fontWeight','bold','fontSize',12);
    set(gca,'fontSize',12,'fontWeight','bold')
    set(findall(gcf,'type','text'),'FontSize',25,'fontWeight','bold'),grid
    
    Pd_ave_save(q) = Pd_ave(1);
end

clear cyc_data

Pcycle_max_ave = max(Pd_ave_save)

