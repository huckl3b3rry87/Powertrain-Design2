%% plot the feasible results

[ROW, COLUMN] = find(FAIL == 0);
L1 = length(ROW);
L2 = length(COLUMN);
S1 = [L1 L2];
S2 = max(S1);
 
MPG_opt = zeros(size(MPG));
 
for num = 1:S2
    MPG_opt(ROW(num),COLUMN(num)) = MPG(ROW(num),COLUMN(num));
end
[A,B,V_min] = find(min(min(MPG(MPG_opt~=0))));
[A,B,V_max] = find(max(max(MPG_opt)));
Vf = [linspace(V_min,V_max,9)];
[FD_opt, MG_opt] = find(V_max == MPG_opt) % (Row, Coulumn) - (y, x)
figure(23);clf
[C,h] = contourf(G_range, FD_range, MPG_opt,Vf);
hold on;
plot(G_range(MG_opt),FD_range(FD_opt),'g.','markersize',40)
% clabel(C,h, 'fontsize', 8);
xlabel('Motor Gear Ratio');
ylabel('Final Drive Ratio');
set(gca,'FontSize',20,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',25,'fontWeight','bold')
legend('MPG','Optimal Gear Ratios')
title(num2str(cyc_name),'fontsize', 30,'fontWeight','bold')
colorbar('YTickLabel',{num2str(Vf(1)),num2str(Vf(2)),num2str(Vf(3)),num2str(Vf(4)),num2str(Vf(5)),num2str(Vf(6))...
    num2str(Vf(7)),num2str(Vf(8)),num2str(Vf(9))},'FontSize',17,'fontWeight','bold')
