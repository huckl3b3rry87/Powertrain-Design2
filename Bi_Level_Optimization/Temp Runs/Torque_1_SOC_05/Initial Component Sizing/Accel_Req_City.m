function [ V_0, V_f,  Acc_Final] = Accel_Req_City( S, dt_2, graph )

dt = 1;
mph_mps = 1/2.237;
C = [[1 1 0];[1 0 1];[0 1 1];[1 0 0];[0 1 0];[0 0 1]]; % Cell array of colros.

for n = 1:3  % Careful with this number for the graph
    cd('Drive_Cycle')
    if n == 1
        %  load CYC_INDIA_URBAN_SAMPLE; cyc_name1 = 'INDIA URBAN';
        % elseif n == 2
        load CYC_MANHATTAN;cyc_name1 = 'MANHATTAN';
    elseif n == 2
        % load CYC_US06; cyc_name3 = 'US06';
        % elseif n == 4
        load CYC_NurembergR36;cyc_name2 = 'NurembergR36';
        %load CONST_65; cyc_name = 'CONST_65';
        %load CONST_45; cyc_name = 'CONST_45';
    elseif n == 3
        load CYC_'NYCC'; cyc_name3 = 'NYCC';
    end
    cd ..
    
    v = cyc_mph(:,2)*mph_mps;  % Cycle Speed in (m/s)
    time_cyc = length(v);
    cyc_time = 1:time_cyc;
    
    % Need to define the a(t) of the profile
    for i = 2:size(v)   % Backwards difference omits first point
        a(i) = (v(i) - v(i-1))/dt;
    end
    
    V_max = 40*mph_mps;
    
    V_test = linspace(0,V_max,S);
    
    % Assume that acceleration is applied at the highest velocity in the interval and increase it by 5% for conservancy
    for i = 1:S-1
        Index = find(V_test(i) < v & v < V_test(i + 1));
        if ~isempty(Index)
            a_req(i) = max(a(Index));
        else
            a_req(i) = 0;
        end
    end
    

    if n~=2  % Exclude Nuremberg
        if n>2      
            a_save(n-1,:) = a_req;
        else
            a_save(n,:) = a_req;
        end
    end
    
    if graph == 1
        figure(1);
        plot(V_test(2:end)/mph_mps,a_req,'color',C(n,:),'linewidth',5)
        hold on
    end
    
end
r =1;
max_a = max(a_save,[],1)*1.05;
for z = 1:length(max_a)
    if max_a(z+1) < max_a(z)
        stop = z;
        break;
    end
end
index = zeros(size(max_a));
for i = (S-1):-1:stop  % Starts decreasing at 3
    if any(max_a(i) < max_a(i+1:S-1))  % Find the next maximum headed backwards and interpolate - except if there is no next maximum - like at the start
        index(i) = 1;
    end
    r = r+1;
end

i = find(index==0);
for g = 1:length(i)
    a_plot(g) = max_a(i(g));
    v_plot(g) = V_test(i(g)+1);
end

Acc_Final = interp1(v_plot,a_plot,V_test(2:end));

if graph == 1
    plot(V_test(2:end)/mph_mps,Acc_Final,'k','linewidth',5)
    legend({cyc_name1,cyc_name2,cyc_name3},'Required Acceleration')
    xlabel('Velocity (MPH)')
    ylabel('Acceleration (m/s^2)')
    set(gca,'fontSize',12,'fontWeight','bold')
    set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','bold'),grid
    hold off
end

V_Initial =  V_test(2:end);
V_Final = (V_Initial + Acc_Final*dt_2);

V_0 = V_Initial/mph_mps;
V_f = V_Final/mph_mps;

end

