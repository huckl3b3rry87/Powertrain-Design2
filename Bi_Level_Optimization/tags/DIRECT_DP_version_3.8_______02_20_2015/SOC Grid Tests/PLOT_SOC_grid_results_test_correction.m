clear;clf;
C = [[1 0 1];[0 1 1];[1 0 0];[0 1 0];[0 0 1]]; % Cell array of colros.
h = 14;
r=6;
for i = 1:5
    
    if i == 1
        load SOC_grid_005_interp;
    elseif i == 2
        load SOC_grid_01;
    elseif i == 3
        load SOC_grid_005;
    elseif i == 4
        load SOC_grid_001;
    else
        load SOC_grid_0005;
    end
    
    
    figure(1)
    subplot(5,1,1)
    plot(result.a1(1,:),result.mpg(1,:),'-ko',...
        'LineWidth',0.5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',C(i,:),...
        'MarkerSize',3)
    grid on
    ylabel('MPG')
    xlabel('NOx weights')
    legend('0.1', '0.01', '0.005','0.001','0.0005');
    hold on
    
    subplot(5,1,2)
    plot(result.a2(2,:),result.mpg(2,:),'-ko',...
        'LineWidth',0.5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',C(i,:),...
        'MarkerSize',3)
    grid on
    ylabel('MPG')
    xlabel('CO weights')
    legend('0.1', '0.01', '0.005','0.001','0.0005');
    hold on
    
    subplot(5,1,3)
    plot(result.a3(3,:),result.mpg(3,:),'-ko',...
        'LineWidth',0.5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',C(i,:),...
        'MarkerSize',3)
    grid on
    ylabel('MPG')
    xlabel('HC weights')
    legend('0.1', '0.01', '0.005','0.001','0.0005');
    hold on
    
    subplot(5,1,4)
    plot(result.shift(4,:),result.mpg(4,:),'-ko',...
        'LineWidth',0.5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',C(i,:),...
        'MarkerSize',3)
    grid on
    ylabel('MPG')
    xlabel('SHIFT weights')
    legend('0.1', '0.01', '0.005','0.001','0.0005');
    hold on
    
    subplot(5,1,5)
    plot(result.eng(5,:),result.mpg(5,:),'-ko',...
        'LineWidth',0.5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',C(i,:),...
        'MarkerSize',3)
    grid on
    ylabel('MPG')
    xlabel('ENG weights')
    legend('0.1', '0.01', '0.005','0.001','0.0005');
    hold on
end

hold off
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','bold')

%%
for i = 1:5
    
    if i == 1
        load SOC_grid_1;
    elseif i == 2
        load SOC_grid_01;
    elseif i == 3
        load SOC_grid_005;
    elseif i == 4
        load SOC_grid_001;
    else
        load SOC_grid_0005;
    end
    
    
    figure(2)
    plot(result.a1(1,:),(result.mpg(1,:)+ result.dSOC(4,:)*48.5),'-ko',...
        'LineWidth',0.5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',C(i,:),...
        'MarkerSize',r)
    grid on
    ylabel('MPG')
    xlabel('NOx weights')
    legend('0.1', '0.01', '0.005','0.001','0.0005');
    set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','bold')
    hold on
    
   figure(3)
    plot(result.a2(2,:),result.mpg(2,:),'-ko',...
        'LineWidth',0.5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',C(i,:),...
        'MarkerSize',r)
    grid on
    ylabel('MPG')
    xlabel('CO weights')
    legend('0.1', '0.01', '0.005','0.001','0.0005');
    set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','bold')
    hold on
    
   figure(4)
    plot(result.a3(3,:),result.mpg(3,:),'-ko',...
        'LineWidth',0.5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',C(i,:),...
        'MarkerSize',r)
    grid on
    ylabel('MPG')
    xlabel('HC weights')
    legend('0.1', '0.01', '0.005','0.001','0.0005');
    set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','bold')
    hold on
    
   figure(5)
    plot(result.shift(4,:),(result.mpg(4,:) + result.dSOC(4,:)*48.5),'-ko',...
        'LineWidth',0.5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',C(i,:),...
        'MarkerSize',r)
    grid on
    ylabel('MPG')
    xlabel('SHIFT weights')
    legend('0.1', '0.01', '0.005','0.001','0.0005');
    set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','bold')
    hold on
    
   figure(6)
    plot(result.eng(5,:),result.mpg(5,:),'-ko',...
        'LineWidth',0.5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',C(i,:),...
        'MarkerSize',r)
    grid on
    ylabel('MPG')
    xlabel('ENG weights')
    legend('0.1', '0.01', '0.005','0.001','0.0005');
    set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','bold')
    hold on
end

hold off
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','bold')

