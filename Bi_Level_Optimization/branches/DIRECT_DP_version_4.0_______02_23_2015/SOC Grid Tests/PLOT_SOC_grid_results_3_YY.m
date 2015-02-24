clear;
clf;
close all;
C = jet;
h = 14;
r=6;
t = 10;
q = 40;  % For colors
n = 4;
h = 9;
y = 12;
x = 1.25;

start = 12;
stop = 10;
for i = [start];
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
    elseif i == 11
        load SOC_grid_005_interp_full_cyc;
        c = '0.005 line. interp FULL';
    else
        load SOC_grid_005_interp_full_cyc_smooth;
        c = '0.005 line. interp FULL Smoothed Data';
    end
    
    if i == start
        leg = {c};
        a = 1;
    else
        leg = {leg{1:a}, c};
        a = a+1;
    end
    
    figure(1)
    [AX,L1,L2] = plotyy(result.a1(1,:),result.mpg(1,:),result.a1(1,:),result.NOx(1,:));
    set(L1,'marker','o',...
        'LineWidth',0.5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',C(a+q,:),...
        'MarkerSize',r);
    set(L2,'marker','o',...
        'LineWidth',0.5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',C(a,:),...
        'MarkerSize',r);
    set(AX, 'fontsize', y,'fontweight','bold');
    ylabel({'MPG'}, 'parent', AX(1));
    ylabel('NOx (g)', 'parent', AX(2));grid
    xlabel('NOx weights')
    title(c);
    H = legend('MPG', 'NOx'); set(H,'fontsize', h)
    
    figure(2)
    [AX,L1,L2] = plotyy(result.a2(2,:),result.mpg(2,:),result.a2(2,:),result.CO(2,:));
    set(L1,'marker','o',...
        'LineWidth',0.5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',C(a+q,:),...
        'MarkerSize',r);
    set(L2,'marker','o',...
        'LineWidth',0.5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',C(a,:),...
        'MarkerSize',r);
    set(AX, 'fontsize', y,'fontweight','bold');
    ylabel({'MPG'}, 'parent', AX(1));
    ylabel('CO (g)', 'parent', AX(2));grid
    xlabel('CO weights')
    title(c);
    H = legend('MPG', 'CO'); set(H,'fontsize', h)
    
    figure(3)
    [AX,L1,L2] = plotyy(result.a3(3,:),result.mpg(3,:),result.a3(3,:),result.HC(3,:));
    set(L1,'marker','o',...
        'LineWidth',0.5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',C(a+q,:),...
        'MarkerSize',r);
    set(L2,'marker','o',...
        'LineWidth',0.5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',C(a,:),...
        'MarkerSize',r);
    set(AX, 'fontsize', y,'fontweight','bold');
    ylabel({'MPG'}, 'parent', AX(1));
    ylabel('HC (g)', 'parent', AX(2));grid
    xlabel('HC weights')
    title(c);
    H = legend('MPG', 'HC'); set(H,'fontsize', h)
    
%         figure(4)
%     [AX,L1,L2] = plotyy(result.shift(4,:),result.mpg(4,:),result.shift(4,:),result.(2,:));
%     set(L1,'marker','o',...
%         'LineWidth',0.5,...
%         'MarkerEdgeColor','k',...
%         'MarkerFaceColor',C(a+q,:),...
%         'MarkerSize',r);
%     set(L2,'marker','o',...
%         'LineWidth',0.5,...
%         'MarkerEdgeColor','k',...
%         'MarkerFaceColor',C(a,:),...
%         'MarkerSize',r);
%     set(AX, 'fontsize', y,'fontweight','bold');
%     ylabel({'MPG'}, 'parent', AX(1));
%     ylabel('CO (g)', 'parent', AX(2));grid
%     xlabel('CO weights')
%     title(c);
%     H = legend('MPG', 'CO'); set(H,'fontsize', h)
    %     figure(3)
    %     plot(result.a2(2,:),result.mpg(2,:),'-ko',...
    %         'LineWidth',0.5,...
    %         'MarkerEdgeColor','k',...
    %         'MarkerFaceColor',C(a+(a-1)*q,:),...
    %         'MarkerSize',r)
    %     grid on
    %     ylabel('MPG')
    %     xlabel('CO weights')
    %     legend(leg{1:a})
    %     set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','bold')
    %     hold on
    %
    %     figure(4)
    %     plot(result.a3(3,:),result.mpg(3,:),'-ko',...
    %         'LineWidth',0.5,...
    %         'MarkerEdgeColor','k',...
    %         'MarkerFaceColor',C(a+(a-1)*q,:),...
    %         'MarkerSize',r)
    %     grid on
    %     ylabel('MPG')
    %     xlabel('HC weights')
    %     legend(leg{1:a})
    %     set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','bold')
    %     hold on
    %
    %     figure(5)
    %     plot(result.shift(4,:),result.mpg(4,:),'-ko',...
    %         'LineWidth',0.5,...
    %         'MarkerEdgeColor','k',...
    %         'MarkerFaceColor',C(a+(a-1)*q,:),...
    %         'MarkerSize',r)
    %     grid on
    %     ylabel('MPG')
    %     xlabel('SHIFT weights')
    %     legend(leg{1:a})
    %     set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','bold')
    %     hold on
    %
    %     figure(6)
    %     plot(result.eng(5,:),result.mpg(5,:),'-ko',...
    %         'LineWidth',0.5,...
    %         'MarkerEdgeColor','k',...
    %         'MarkerFaceColor',C(a+(a-1)*q,:),...
    %         'MarkerSize',r)
    %     grid on
    %     ylabel('MPG')
    %     xlabel('ENG weights')
    %     legend(leg{1:a})
    %     set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','bold')
    %     hold on
end

hold off
set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','bold')

