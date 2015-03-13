clear;
clf;
close all;
C = jet;
h = 14;
r=6;
t = 10;
q = 20;  % For colors
start = 3;
stop = 10;
for i = [start, 6,8,9];
    if i == 1
        load SOC_grid_1;
        c = '0.1';
    elseif i == 2
        load SOC_grid_01;
        c = '0.01';
    elseif i == 3
        load SOC_grid_005;
        c = '0.005';
    elseif i == 4
        load SOC_grid_001;
        c = '0.001';
    elseif i == 5
        load SOC_grid_0005;
        c = '0.0005';
    elseif i == 6
        load SOC_grid_005_interp;
        c = '0.005 lin. interp';
    elseif i == 7
        load SOC_grid_001_interp;
        c = '0.001 lin. interp.';
    elseif i == 8
        load SOC_grid_0005_interp;
        c = '0.0005 lin. interp.';
    elseif i == 9
        load SOC_grid_005_interp_cubic;
        c = '0.005 cubic interp.';
    elseif i == 10
        load SOC_grid_0005_interp_cubic;
        c = '0.0005 cubic. interp.';
    elseif i == 11;
        load SOC_grid_005_interp_full_cyc;
        c = '0.005 line. interp FULL';
    elseif i == 12;
        load SOC_grid_005_interp_lower_weights;
        c = '0.005 line. interp Lower weights';
    elseif i == 13;
        load SOC_grid_005_interp_full_cyc_smooth;
        c = '0.005 line. interp Full CYC Smoothed Data';
    else
        load SOC_grid_001_interp_full_cyc_smooth;
        c = '0.001 line. interp Full CYC Smoothed Data';
    end
    
    if i == start
        leg = {c};
        a = 1;
    else
        leg = {leg{1:a}, c};
        a = a+1;
    end
    
    figure(2)
    plot(result.a1(1,:),result.mpg(1,:),'-ko',...
        'LineWidth',0.5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',C(a+(a-1)*q,:),...
        'MarkerSize',r)
    grid on
    ylabel('MPG')
    xlabel('NOx weights')
    legend(leg{1:a})
    set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','bold')
    hold on
    
    figure(3)
    plot(result.a2(2,:),result.mpg(2,:),'-ko',...
        'LineWidth',0.5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',C(a+(a-1)*q,:),...
        'MarkerSize',r)
    grid on
    ylabel('MPG')
    xlabel('CO weights')
    legend(leg{1:a})
    set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','bold')
    hold on
    
    figure(4)
    plot(result.a3(3,:),result.mpg(3,:),'-ko',...
        'LineWidth',0.5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',C(a+(a-1)*q,:),...
        'MarkerSize',r)
    grid on
    ylabel('MPG')
    xlabel('HC weights')
    legend(leg{1:a})
    set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','bold')
    hold on
    
    figure(5)
    plot(result.shift(4,:),result.mpg(4,:),'-ko',...
        'LineWidth',0.5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',C(a+(a-1)*q,:),...
        'MarkerSize',r)
    grid on
    ylabel('MPG')
    xlabel('SHIFT weights')
    legend(leg{1:a})
    set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','bold')
    hold on
    
    figure(6)
    plot(result.eng(5,:),result.mpg(5,:),'-ko',...
        'LineWidth',0.5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',C(a+(a-1)*q,:),...
        'MarkerSize',r)
    grid on
    ylabel('MPG')
    xlabel('ENG weights')
    legend(leg{1:a})
    set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','bold')
    hold on
end

hold off
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','bold')

