clear
clc
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%-----------------Define the Run Type-------------------------------------%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
RUN_TYPE.sim = 0;  % RUN_TYPE = 1 - for DIRECT     &    RUN_TYPE = 0 - for DP only
RUN_TYPE.emiss = 1; % RUN_TYPE.emiss = 1 - has emissions  &   RUN_TYPE.emiss = 0 - NO emissions
RUN_TYPE.plot = 1;  % RUN_TYPE.plot = 1 - plots on  &   RUN_TYPE.plot = 0 - plots off
RUN_TYPE.soc_size = 0.05;
RUN_TYPE.trq_size = 1;  % Nm

%% eval(['save(''',''SOC_grid_005_trq_20_NM_HWFET',''',','''result'');'])
%  eval(['save(''','SOC_grid_05_trq_1_NM_HWFET',''',','''result'');'])

% Code Parameters
beta = 20;
alpha = 20;
h = 20;
if RUN_TYPE.emiss == 1
    for k = 1:5
        for u = 1:alpha
            if k ==1;
                a1=(u-1)/beta;
                a2=0;
                a3=0;
                shift = 0;
                eng = 0;
            elseif k==2
                a1=0;
                a2 = (u-1)/beta;
                a3=0;
                shift = 0;
                eng = 0;
            elseif k ==3
                a1=0;
                a2 = 0;
                a3=(u-1)/beta;
                shift = 0;
                eng = 0;
            elseif k ==4
                a1=0;
                a2=0;
                a3=0;
                shift = (u-1)/beta;
                eng= 0;
            else
                a1=0;
                a2=0;
                a3=0;
                shift = 0;
                eng = (u-1)/beta;
            end
            
            [ FAIL, MPG, emission, delta_SOC, sim ] = find_weights(a1,a2,a3,shift,eng,RUN_TYPE);
            
            if FAIL.final == 1;
                result.mpg(k,u) = NaN;
                result.NOx(k,u) = NaN;
                result.HC(k,u)= NaN;
                result.CO(k,u) = NaN;
                result.EE(k,u) = NaN;
                result.SE(k,u) = NaN;
            else
                result.mpg(k,u) = MPG;
                result.NOx(k,u) = emission.NOx;
                result.HC(k,u)= emission.HC;
                result.CO(k,u) = emission.CO;
                result.EE(k,u) = sim.EE;
                result.SE(k,u) = sim.SE;
            end
            result.dSOC(k,u) = delta_SOC;
            result.a1(k,u) = a1;
            result.a2(k,u) = a2;
            result.a3(k,u) = a3;
            result.shift(k,u) = shift;
            result.eng(k,u) = eng;
            result.feasible(k,u) = FAIL; % If it is zero then it is OK
            plot(u,k,'r.','markersize',25)
            ylabel('k total ')
            xlabel('u - 20 per'),grid
            hold on
        end
    end
    % Total = [result.mpg; result.NOx; result.HC; result.CO; result.dSOC]
    clf;
    figure(1)
    subplot(5,1,1)
    plot(result.a1(1,:),result.mpg(1,:),'r.','markersize',h)
    hold on
    plot(result.a1(1,:),result.mpg(1,:))
    ylabel('MPG')
    xlabel('NOx weights'),grid
    
    subplot(5,1,2)
    plot(result.a2(2,:),result.mpg(2,:),'r.','markersize',h)
    hold on
    plot(result.a2(2,:),result.mpg(2,:))
    ylabel('MPG')
    xlabel('CO weights'),grid
    
    subplot(5,1,3)
    plot(result.a3(3,:),result.mpg(3,:),'r.','markersize',h)
    hold on
    plot(result.a3(3,:),result.mpg(3,:))
    ylabel('MPG')
    xlabel('HC weights'),grid
    
    subplot(5,1,4)
    plot(result.shift(4,:),result.mpg(4,:),'r.','markersize',h)
    hold on
    plot(result.shift(4,:),result.mpg(4,:))
    ylabel('MPG')
    xlabel('SHIFT weights'),grid
    
    subplot(5,1,5)
    plot(result.eng(5,:),result.mpg(5,:),'r.','markersize',h)
    hold on
    plot(result.eng(5,:),result.mpg(5,:))
    ylabel('MPG')
    xlabel('ENG weights'),grid
    
else  %% No emissions
    for k = 1:2
        for u = 1:alpha
            if k == 1
                a1=0;
                a2=0;
                a3=0;
                shift = (u-1)/beta;
                eng= 0;
            else
                a1=0;
                a2=0;
                a3=0;
                shift = 0;
                eng = (u-1)/beta;
            end
            
            [ FAIL, MPG, emission, delta_SOC ] = find_weights(a1,a2,a3,shift,eng)
            
            if FAIL == 1;
                result.mpg(k,u) = NaN;
            else
                result.mpg(k,u) = MPG;
            end
            result.dSOC(k,u) = delta_SOC;
            result.a1(k,u) = a1;
            result.a2(k,u) = a2;
            result.a3(k,u) = a3;
            result.shift(k,u) = shift;
            result.eng(k,u) = eng;
        end
    end
    figure(1)
    subplot(2,1,1)
    plot(result.shift(1,:),result.mpg(1,:),'r.','markersize',h)
    hold on
    plot(result.shift(1,:),result.mpg(1,:))
    ylabel('MPG')
    xlabel('SHIFT weights'),grid
    
    subplot(2,1,2)
    plot(result.eng(2,:),result.mpg(2,:),'r.','markersize',h)
    hold on
    plot(result.eng(2,:),result.mpg(2,:))
    ylabel('MPG')
    xlabel('ENG weights'),grid
end
 eval(['save(''','SOC_grid_05_trq_1_NM_HWFET',''',','''result'');'])
