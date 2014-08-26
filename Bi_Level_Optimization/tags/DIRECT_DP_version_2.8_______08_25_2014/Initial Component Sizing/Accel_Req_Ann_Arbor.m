function [ V_0, V_f,  Acc_Final] = Accel_Req_Ann_Arbor( S, dt_2, graph )

dt = 1;
mph_mps = 1/2.237;
C = [[1 1 0];[1 0 1];[0 1 1];[0.752807, 0.752907, 0.752907];[0 1 0];[0 0 1]; [75,0,130]/130;[1, 0.4, 0]; [51,102,0]/102; [290,0.4,0.6]/290]; % Cell array of colros.
g = 1;
for n = 1:10  % Careful with this number for the graph
    cd('Drive_Cycle')
    if n == 1
        load D1;cyc_name1 = 'Driver 1';
    elseif n == 2
        load D2;cyc_name2 = 'Driver 2';
    elseif n == 3
        load D3; cyc_name3 = 'Driver 3';
    elseif n == 4
        load D4;cyc_name4 = 'Driver 4';
    elseif n == 5
        load D5; cyc_name5 = 'Driver 5';
    elseif n == 6
        load D6;cyc_name6 = 'Driver 6';
    elseif n == 7
        load D7; cyc_name7 = 'Driver 7';
    elseif n == 8
        load D8;cyc_name8 = 'Driver 8';
    elseif n == 9
        load D9; cyc_name9 = 'Driver 9';
    else % n == 10
        load D10; cyc_name10 = 'Driver 10';
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
    
    if n~=1 && n~=2 && n~= 9 && n~=6 && n~=7 && n~=3  
        a_save(g,:) = a_req;
        g = g + 1;
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
    legend({cyc_name1,cyc_name2,cyc_name3,cyc_name4,cyc_name5,cyc_name6,cyc_name7,cyc_name8,cyc_name9,cyc_name10},'Required Acceleration')
    xlabel('Velocity (mph)')
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

